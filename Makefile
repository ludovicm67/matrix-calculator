tp: tp.c
	gcc tp.c -lm -o tp

.PHONY: clean
clean:
	rm -f tp
