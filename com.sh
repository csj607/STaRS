mpiifort -ipo -O3 -no-prec-div -fp-model fast=2 -cpp -DMPI -c random_mt.f90
mpiifort -ipo -O3 -no-prec-div -fp-model fast=2 -cpp -DMPI -c RT*.f90
mpiifort -ipo -O3 -no-prec-div -fp-model fast=2 -cpp -DMPI -c memory_mod.f90
mpiifort -ipo -O3 -no-prec-div -fp-model fast=2 -cpp -DMPI main.f90 *.o
