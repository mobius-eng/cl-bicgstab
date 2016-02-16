(asdf:defsystem #:cl-bicgstab
  :description "BiCGStab method implementation"
  :depends-on (:cl-linear-algebra :cl-iterative :cl-numerics-utils)
  :serial t
  :components ((:file "package")
               (:file "bicgstab")))
