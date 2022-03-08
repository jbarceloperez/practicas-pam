all:
	gcc -o mulmatcuablo MulMatCuaBlo.c -O3 -lm -fopenmp
	gcc -o secmulmatblo SecMulMatBlo.c -O3 -lm -fopenmp
	gcc -o secmulmatcuablo SecMulMatCuaBlo.c -O3 -lm -fopenmp
	gcc -o mulmatblo MulMatBlo.c -O3 -lm -fopenmp

mulmatcuablo: MulMatCuaBlo.c
	gcc -o mulmatcuablo MulMatCuaBlo.c -O3 -lm -fopenmp	

secmulmatblo: SecMulMatBlo.c
	gcc -o secmulmatblo SecMulMatBlo.c -O3 -lm -fopenmp	

secmulmatcuablo: SecMulMatCuaBlo.c
	gcc -o secmulmatcuablo SecMulMatCuaBlo.c -O3 -lm -fopenmp
mulmatblo: MulMatBlo.c
	gcc -o mulmatblo MulMatBlo.c -O3 -lm -fopenmp

clear:
	rm mulmatcuablo secmulmatblo secmulmatcuablo mulmatblo
