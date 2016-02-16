(in-package cl-user)

(defpackage :cl-bicgstab
  (:use #:cl #:cl-linear-algebra #:cl-iterative #:cl-numerics-utils #:optima)
  (:export #:bicg-stab #:bicg-stab-size #:bicg-stab-residual #:bicg-stab-solution
           #:bicg-stab-value
           #:bicg-stab-solve))
