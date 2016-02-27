(in-package cl-bicgstab)
;; * BiCGStab method implementation

;; ** Value: state of BiCGStab algorithm
(defclass bicg-stab-value ()
  ((bicg-rho :initform 1.0d0 :accessor bicg-rho)
   (bicg-alpha :initform 1.0d0 :accessor bicg-alpha)
   (bicg-omega :initform 1.0d0 :accessor bicg-omega)
   (bicg-v :initarg :v :reader bicg-v)
   (bicg-p :initarg :p
           :accessor bicg-p
           :documentation "Search direction for next approximation")
   (bicg-r0 :initarg :r0
            :reader bicg-r0
            :documentation "Residual of the initial approximation")
   (bicg-r :initarg :r
           :reader bicg-r
           :documentation "Residual of approximation")
   (bicg-x :initarg :x
           :reader bicg-x
           :documentation "Current approximation")
   (bicg-s :initarg :s :reader bicg-s)
   (bicg-t :initarg :t :reader bicg-t))
  (:documentation "State of BiCGStab computation"))

(defmethod print-object ((obj bicg-stab-value) out)
  (print-unreadable-object (obj out :type nil)
    (format out "BICG-STAB~%x = ~A~%r = ~A~%r0 = ~A~%"
            (bicg-x obj) (bicg-r obj) (bicg-r0 obj))))

(defun bicg-stab-value (size)
  "Make BICG-DATA with all vector sizes SIZE"
  (make-instance 'bicg-stab-value
    :v (make-vector size 'double-float)
    :p (make-vector size 'double-float)
    :r0 (make-vector size 'double-float)
    :r (make-vector size 'double-float)
    :x (make-vector size 'double-float)
    :s (make-vector size 'double-float)
    :t (make-vector size 'double-float)))

(defun bicg-stab-vector-length (bicg-stab-value)
  "Size of BiCGStab problem"
  (length (bicg-x bicg-stab-value)))

(defun residual! (A x b r)
  "Calculate the residual r = b - A(x)"
  (declare (optimize (speed 3) (safety 1) (debug 1)))
  (declare (type (simple-array double-float *) r)
           (type (simple-array * *) x b)
           (type function A))
  (funcall A x r)
  (dotimes (i (length r))
    (setf (aref r i) (- (aref b i) (aref r i)))))

(defun init-value! (value A x0 b)
  "Initialize value with approximation X0
Residual (BICG-R) is initialized as b - A*x0"
  (declare (optimize (speed 3) (safety 1) (debug 0)))
  (declare (type bicg-stab-value value))
  (setf (bicg-alpha value) 1.0d0)
  (setf (bicg-rho value) 1.0d0)
  (setf (bicg-omega value) 1.0d0)
  (residual! A x0 b (bicg-r0 value))
  (zeros! (bicg-p value))
  (zeros! (bicg-v value))
  (copy-vector-contents (bicg-r0 value) (bicg-r value))
  (copy-vector-contents x0 (bicg-x value)))

;; ** BiCGStab computation steps
(defun new-rho ()
  " Calculates rho = (dot r r0)."
  (alter-value
   (lambda (val) (cons (dot (bicg-r val) (bicg-r0 val)) val))))

(defun restart-if-rho-is-zero (tolerance)
  "Cheks if new-rho is too close to zero and restarts computation
from current approximation. I.e. r0 is set to current residual."
  (alter-value
   (lambda (x)
     (declare (optimize (speed 3) (debug 0) (safety 1)))
     (destructuring-bind (new-rho . value) x
       (if (num= 0d0 new-rho tolerance)
           (progn
             (copy-vector-contents (bicg-r value) (bicg-r0 value))
             (zeros! (bicg-v value))
             (zeros! (bicg-p value))
             (setf (bicg-rho value) 1d0)
             (setf (bicg-alpha value) 1d0)
             (setf (bicg-omega value) 1d0)
             (cons (dot (bicg-r value) (bicg-r0 value)) value))
           x)))))

(defun new-direction (A)
  "Get new direction

    beta = rho    / rho  * alpha / omega
              i-1      i                i-1
 
    p   =  r    +  beta (p   - omega    v   )
     i      i-1           i-1       i-1  i-1 

   v  = A p
    i      i

"
  (declare (type function A))
  (alter-value
   (lambda (x)
     (declare (optimize (speed 3) (safety 1) (debug 0)))
     (destructuring-bind (new-rho . v) x
       (declare (type double-float new-rho))
       (let ((beta (* (/ new-rho (the double-float (bicg-rho v)))
                      (/ (the double-float (bicg-alpha v))
                         (the double-float (bicg-omega v))))))
         (declare (type double-float beta))
         (lincomb beta (bicg-p v)
                  :vectors (list (bicg-r v) (bicg-v v))
                  :multipliers (list 1d0 (- (* beta (the double-float (bicg-omega v))))))
         (funcall A (bicg-p v) (bicg-v v))
         (setf (bicg-rho v) new-rho)
         v)))))

(defun r0*v ()
  "Calculate dot-product (r0, v)"
  (alter-value
   (lambda (value)
     (cons (dot (bicg-r0 value) (bicg-v value) nil) value))))

;; Should never fail here, but just in case
(defun r0*v-is-not-zero (tolerance)
  "Fail if r0*v is zero"
  (failed-value
   (lambda (x)
     (let ((r0*v (car x)))
       (num= 0d0 r0*v (expt tolerance 2))))
   (lambda (x)
     (list "r0*v is zero" (cdr x)))))

(defun new-alpha-s ()
  " Calculate new alpha:

               rho
                  i
   alpha = -----------
             r0 * v
                   i

and

  s = r    - alpha * v
       i-1            i
"
  (alter-value
   (lambda (x)
     (destructuring-bind (r0*v . v) x
       (setf (bicg-alpha v) (/ (bicg-rho v) r0*v))
       (sum-vectors
        (bicg-s v)
        :vectors (list (bicg-r v) (bicg-v v))
        :multipliers (list 1d0 (- (bicg-alpha v))))
       v))))

(defun finish-if-s-is-small (ewt)
  "If rms-norm of s is small - finish"
  (finished-value
   (lambda (bicg-value)
     (let ((s-norm (rms-norm (bicg-s bicg-value) (ewt-vector ewt (bicg-x bicg-value)))))
       (<= s-norm 1d0)))
   (lambda (value)
     (add-vectors (bicg-x value)
                  :vectors (list (bicg-p value))
                  :multipliers (list (bicg-alpha value)))
     (copy-vector-contents (bicg-s value) (bicg-r value))
     value)))

(defun new-t-omega-x-r (A)
  "Update:

    t = A s

               (s, t)  
    omega  = -----------
         i     (t, t)


    x  = x    + alpha p  + omega  s
     i    i-1          i        i

    r  = s - omega  t
     i            i

"
  (alter-value
   (lambda (value)
     (funcall A (bicg-s value) (bicg-t value))
     (setf (bicg-omega value)
           (/ (dot (bicg-s value) (bicg-t value) nil)
              (square-vector (bicg-t value) nil)))
     (add-vectors
      (bicg-x value)
      :vectors (list (bicg-p value) (bicg-s value))
      :multipliers (list (bicg-alpha value) (bicg-omega value)))
     (sum-vectors
      (bicg-r value)
      :vectors (list (bicg-s value) (bicg-t value))
      :multipliers (list 1d0 (- (bicg-omega value))))
     value)))

(defun finish-if-residual-is-small (ewt)
  "Final finish check: finish is residual is small"
  (finished-value
   (lambda (value)
     (let ((r-norm (rms-norm (bicg-r value) (ewt-vector ewt (bicg-x value)))))
       (<= r-norm 1d0)))))

;; *** Computation parameters
(defvar *bicg-stab-tolerance* (sqrt double-float-epsilon))
(defvar *bicg-stab-max-iter-coeff* 10)

(defun bicg-stab-solve (value A b x0 ewt &rest other-controls)
  "Solve linear set of equations A*x=b using BiCGStab method
  VALUE is a compatible in size instance of BICG-STAB-VALUE
  A is a function representing matrix multiplication or linear operator
    (funcall A x b)
  B is a right hand side vector
  X0 is initial approximation of the solution
  EWT is error-weight object containing relative and absolute tolerances
  OTHER-CONTROLS extra controls (e.g. log-value) for computation"
  (check-vector-lengths b x0)
  (assert (= (length b) (bicg-stab-vector-length value))
          ()
          'vector-length-mismatch
          :length1 (length b)
          :length2 (bicg-stab-vector-length value))
  (init-value! value A x0 b)
  (let ((n (length b)))
    (let ((computation (combine-controls
                        (new-rho) (restart-if-rho-is-zero *bicg-stab-tolerance*)
                        (new-direction A)
                        (r0*v) (r0*v-is-not-zero *bicg-stab-tolerance*)
                        (new-alpha-s) (finish-if-s-is-small ewt)
                        (new-t-omega-x-r A)
                        (finish-if-residual-is-small ewt)
                        (limit-iterations (* n *bicg-stab-max-iter-coeff*)
                                          (lambda (n) (list "exceeded max iterations" n))))))
      (iterate (iterator:continue value)
               (apply #'combine-controls computation other-controls)
               (combine-controls
                (finish-if-residual-is-small ewt)
                (finished-value
                 (lambda (value)
                   (declare (ignore value))
                   (<= (rms-norm b (ewt-vector ewt x0)) 1d0))
                 (lambda (value)
                   (zeros! (bicg-x value))
                   (zeros! (bicg-r value))
                   value)))))))

(defclass bicg-stab ()
  ((bicg-value
    :initarg :value
    :reader bicg-value
    :documentation "Value part of BiCGStab method")
   (bicg-other-controls
    :initarg :other-controls
    :reader bicg-other-controls
    :documentation "Other controls for solver")
   (bicg-ewt
    :initarg :ewt
    :reader bicg-ewt
    :documentation "Error weight object for BiCGStab method"))
  (:documentation
   "BiCGStab method"))

(defun bicg-stab (size
                  &optional (rtol *bicg-stab-tolerance*) (atol *bicg-stab-tolerance*)
                  &rest other-controls)
  "Make BiCGSTab method for a problem of size SIZE."
  (make-instance 'bicg-stab
    :value (bicg-stab-value size)
    :other-controls other-controls
    :ewt (ewt (if (vectorp rtol)
                  rtol
                  (make-vector size 'double-float rtol))
              (if (vectorp atol)
                  atol
                  (make-vector size 'double-float atol)))))

(defun bicg-stab-size (bicg-stab)
  "Method's problem size (vectors lengths)"
  (bicg-stab-vector-length (bicg-value bicg-stab)))

(defun bicg-stab-residual (value)
  "Get residual (error) from BICG-STAB-VALUE"
  (declare (type bicg-stab-value value))
  (bicg-r value))

(defun bicg-stab-solution (value)
  "Get solution from BICG-STAB-VALUE"
  (declare (type bicg-stab-value value))
  (bicg-x value))

(defmethod solve-linear ((method bicg-stab) a b x)
  (match (apply #'bicg-stab-solve
                (bicg-value method)
                a
                b
                x
                (bicg-ewt method)
                (bicg-other-controls method))
    ((iterator:iterator :status :finished :value value)
     (copy-vector-contents (bicg-x value) x)
     (values x t (rms-norm (bicg-r value) (ewt-vector (bicg-ewt method) (bicg-x value)))))
    ((iterator:iterator :status :failed :value x)
     (let ((msg (first x))
           (value (second x)))
       (values x nil (bicg-r value) msg)))))


