## Test environments

* local R installation, R 4.5.2
* check_win_devel()

## R CMD check results

0 errors | 1 warnings | 1 notes
W  checking whether package ‘ridgetorus’ can be installed (6.6s)
   Found the following significant warnings:
     /Library/Frameworks/R.framework/Resources/include/R_ext/Boolean.h:62:36: warning: unknown warning group '-Wfixed-enum-extension', ignored [-Wunknown-warning-option]
N  checking for future file timestamps (1m 0.8s)
   unable to verify current time

## Comments

Warning and note are benign and do not show up in the check_win_devel() checks