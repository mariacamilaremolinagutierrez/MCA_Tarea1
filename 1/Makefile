all: plots.py
	python plots.py

plots.py: shocktube.c
	gcc shocktube.c -lm -o shocktube.x
	./shocktube.x
