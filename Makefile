task:
	echo "Labwork 1"
	mpicc main.c -o build/task.out
	mpirun -np 4 --oversubscribe build/task.out
clean:
	rm build/task.out
	rm output.csv
