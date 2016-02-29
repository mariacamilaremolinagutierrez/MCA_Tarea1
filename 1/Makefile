all: plots.py
	python plots.py

plots.py: final_step.dat
	./shocktube.x

final_step.dat: shocktube.c
	gcc shocktube.c -lm -o shocktube.x

clean:
	rm -r Graficas/
	rm shocktube.x
	rm final_step.dat
