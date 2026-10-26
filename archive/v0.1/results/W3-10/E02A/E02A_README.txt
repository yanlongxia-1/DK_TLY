E02A compares MVMD convergence on the SAME E02 validation signal.
K and alpha remain fixed, no retuning.
AC_reconstruction_relative = norm(x-mean(x)-sum(IMFs))/norm(x-mean(x)).
This matters more than residual normalized by x including large sensor offsets.
E02A cannot establish true SNR or real angular-rate accuracy.
Do not choose optimum based on validation SNR or phase results; use diagnostic only.
Send E02A_convergence.csv, E02A_modes.csv, E02A_diagnostic.png.
