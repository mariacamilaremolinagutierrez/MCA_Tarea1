all: Das_Executor.x Outputs/int_RK4_0.353553.dat Outputs/int_simpletico_0.353553.dat Outputs/int_simpletico_0.432500.dat Outputs/int_simpletico_0.425000.dat Outputs/int_simpletico_0.460000.dat Images/Figura_9.pdf Images/Figura_10.pdf Images/Figura_11.pdf Images/Figura_12.pdf Images/Figura_13.pdf Images/Figura_14.pdf Images/Figura_15.pdf Images/Figura_16.pdf Images/Figura_17.pdf

Images/Figura_9.pdf:Outputs/int_RK4_0.353553.dat Outputs/int_simpletico_0.353553.dat graficar.py
	python graficar.py 9 Outputs/int_simpletico_0.353553.dat Outputs/int_RK4_0.353553.dat
Images/Figura_10.pdf:Outputs/int_RK4_0.353553.dat Outputs/int_simpletico_0.353553.dat graficar.py
	python graficar.py 10 Outputs/int_simpletico_0.353553.dat Outputs/int_RK4_0.353553.dat
Images/Figura_11.pdf:Outputs/int_RK4_0.353553.dat Outputs/int_simpletico_0.353553.dat graficar.py
	python graficar.py 11 Outputs/int_simpletico_0.353553.dat Outputs/int_RK4_0.353553.dat
Images/Figura_12.pdf: Outputs/int_simpletico_0.432500.dat graficar.py
	python graficar.py 12 Outputs/int_simpletico_0.432500.dat 
Images/Figura_13.pdf: Outputs/int_simpletico_0.432500.dat graficar.py
	python graficar.py 13 Outputs/int_simpletico_0.432500.dat 
Images/Figura_14.pdf: Outputs/int_simpletico_0.432500.dat graficar.py
	python graficar.py 14 Outputs/int_simpletico_0.432500.dat 
Images/Figura_15.pdf: Outputs/int_simpletico_0.425000.dat graficar.py
	python graficar.py 15 Outputs/int_simpletico_0.425000.dat
Images/Figura_16.pdf: Outputs/int_simpletico_0.460000.dat graficar.py
	python graficar.py 16 Outputs/int_simpletico_0.460000.dat
Images/Figura_17.pdf:Outputs/int_RK4_0.353553.dat Outputs/int_simpletico_0.353553.dat graficar.py
	python graficar.py 17 Outputs/int_simpletico_0.353553.dat Outputs/int_RK4_0.353553.dat
Outputs/int_RK4_0.353553.dat:Das_Executor.x
	./Das_Executor.x 2 2800 0.6 0.35355339059 ./ICs/ics_fig_9_10_111.dat
Outputs/int_simpletico_0.353553.dat:Das_Executor.x
	./Das_Executor.x 1 2800 0.6 0.35355339059 ./ICs/ics_fig_9_10_111.dat
Outputs/int_simpletico_0.432500.dat:Das_Executor.x
	./Das_Executor.x 1 2800 0.6 0.4325 ./ICs/ics_fig_12_13_142.dat
Outputs/int_simpletico_0.425000.dat:Das_Executor.x
	./Das_Executor.x 1 2800 0.6 0.425 ./ICs/ics_fig_153.dat
Outputs/int_simpletico_0.460000.dat:Das_Executor.x
	./Das_Executor.x 1 2800 0.6 0.46 ./ICs/ics_164.dat

Das_Executor.x:Das_Punto_2.c
	cc Das_Punto_2.c -lm -o Das_Executor.x

clean:
	rm /Outputs/*.dat Images/*.pdf
