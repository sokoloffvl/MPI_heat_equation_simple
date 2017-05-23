// ConsoleApplication2.cpp:
//
//mpicc -o final final.c -lm
//mpirun -np 2 final
#include <stdio.h> 
#include <stdint.h>
#include <stdlib.h>
#include <math.h> 
#include <mpi.h>
//Основные глобавльные неизменяемые переменные
#define X_0 0
#define X_N 1
#define Y_0 0
#define Y_N 0.5
#define Z_0 0
#define Z_N 2
#define T_0 0
#define T_N 10

#define sigma 0.01

int max_range;

//Шаги по координатам
double h_x;
double h_y;
double h_z;
//тао :)
double tao;
//расчет значений x,y,z в i-м узле сетки
double x_(int i)
{
	return X_0 + h_x * i;
}

double y_(int j)
{
	return Y_0 + h_y * j;
}

double z_(int k)
{
	return Z_0 + h_z * k;
}

//Функция расчета следущего узла
double temperature_t(int i, int j, int k, int t, double arr[][max_range][max_range])
{
	double result = -1;
		if (i == 0)
		{
			result = 500.0;
		}
		else if (i == max_range - 1)
		{
			result = 1000.0;
		}
		else if (j == 0)
		{
			result = 500.0 * (1 + x_(i));
		}
		else if (j == max_range - 1)
		{
			result = 500.0 * (1 + pow(x_(i), 2));
		}
		else if (k == 0)
		{
			result = 500.0;
		}
		else if (k == max_range - 1)
		{
			result = 500.0 * (1 + pow(x_(i), 4));
		}
		else
		{
			if (t == 0) result = 300.0;
			else result = arr[i][j][k] + sigma*tao*(
					((arr[i-1][j][k] + arr[i+1][j][k] - 2*arr[i][j][k])/(pow(h_x,2)))
					+((arr[i][j-1][k] + arr[i][j+1][k] - 2*arr[i][j][k])/(pow(h_y,2)))
					+((arr[i][j][k-1] + arr[i][j][k+1] - 2*arr[i][j][k])/(pow(h_z,2)))
				);
		}
	return result;
}
//вывод матрицы на экран(потом переделаю для вывода в файл)
void print_matrix(double arr[][max_range][max_range], double spentTime, char* fileName)
{	
	int i, j, k;
	FILE * file;
	file = fopen(fileName, "w");
	fprintf(file,"summary spent time:%3f\n", spentTime);
	fprintf(file,"----------------------------\n");
	fprintf(file, "Final result after all iteartions \n");					
	for (i = 0; i < max_range; i++)
	{
		fprintf(file,"\n");
		fprintf(file,"layer %3d \n", i);
		fprintf(file,"\n");
		for (j = 0; j < max_range; j++)
		{
			for (k = 0; k < max_range; k++)
			{
				fprintf(file,"%3f ", arr[i][j][k]);
			}
			fprintf(file,"\n");
		}
	}
	fprintf(file,"\n");
	fclose(file);
}

void print_temperature(double arr[][max_range][max_range], int t, char* fileName)
{	
	int i, j, k;
	FILE * file;
	file = fopen(fileName, "w");
	fprintf(file, "graph of t = %3d\n", t);	
	fprintf(file,"----------------------------\n");				
	for (i = 0; i < max_range; i++)
	{
		fprintf(file,"%3f\n", arr[i][max_range/2][max_range/2]);
	}
	fprintf(file,"\n");
	fclose(file);
}

int main(int argc, char* argv[])
{

	int i, j, k, t;
	int errCode;

	max_range = 50;

	double M_2[max_range][max_range][max_range];
	double M_1[max_range][max_range][max_range];

	double t_start, t_end;

	h_x = (double)fabs(X_N - X_0) / (double)max_range;
	h_y = (double)fabs(Y_N - Y_0) / (double)max_range;
	h_z = (double)fabs(Z_N - Z_0) / (double)max_range;

	tao = (1 / (2*sigma*((1/pow(h_x,2)) + (1/pow(h_y,2)) + (1/pow(h_z,2)))));
	double delta_t = 1.0;//для удобства дебага. Если сделать delta_t = 1.0; То считаться будет быстро, так как меньше узлов будет по времени
	double t_max = ceil((T_N-T_0)/((double)delta_t));

	char* file1 = "my_2.txt";
	char* file2 = "my_4.txt";
	char* file3 = "my_6.txt";
	char* file4 = "my_10.txt";
	char* file5 = "final.txt";

	MPI_Init(&argc, &argv);// инициализация mpi
	int rank, count, work_count;

	MPI_Status status1, status2;
	MPI_Request req1,req2;

	MPI_Request reqs[T_N*max_range*max_range];

	MPI_Comm_size(MPI_COMM_WORLD, &count);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//левая и правая граница ответственности каждого процесса
	int i_left, i_right;
	//число "считающих процессов", нулевой процесс будет дерижировать
	work_count = count - 1;

	if (rank > 0)
	{
		//определяем границы для процесса № rank
		i_left = 0 + (max_range / work_count) * (rank - 1);
		i_right = 0 + (max_range / work_count) * rank;

		if (rank == count - 1 && i_right < max_range)
			i_right = max_range; //поправляем размер слоя для у процесса с последним слоем

		//printf("process # %3d , my left = %3d, my right = %3d \n", rank, i_left, i_right);

		//Расчет матрицы
		for (t = T_0; t <= t_max; t++)
		{
			if (t == 0)
			{
				for (i = i_left; i < i_right; i++)
					for (j = 0; j < max_range; j++)
						for (k = 0; k < max_range; k++)
							M_1[i][j][k] = temperature_t(i, j, k, t, M_2);
				//отправляем направо
				if (rank < count - 1)
				{
					MPI_Isend(&M_1[i_right-1][0][0], max_range*max_range, MPI_DOUBLE, rank+1, (t+1)*(rank*count+(rank+1)), MPI_COMM_WORLD, &req1);
				}
				//отправляем налево
				if (rank > 1)
				{
					MPI_Isend(&M_1[i_left][0][0], max_range*max_range, MPI_DOUBLE, rank-1, (t+1)*(rank*count+(rank-1)), MPI_COMM_WORLD, &req2);
				}
			}
			else
			{
				//получаем слева
				if (rank > 1)
				{
					MPI_Irecv(&M_2[i_left-1][0][0], max_range*max_range, MPI_DOUBLE, rank-1, t*((rank-1)*count+rank), MPI_COMM_WORLD, &req1);
					MPI_Wait(&req1, &status1);
				}
				//получаем справа
				if (rank < count - 1)
				{
					MPI_Irecv(&M_2[i_right][0][0], max_range*max_range, MPI_DOUBLE, rank+1, t*((rank+1)*count+rank), MPI_COMM_WORLD, &req2);
					MPI_Wait(&req2, &status2);
				}
				//считаем в M1
				for (i = i_left; i < i_right; i++)
					for (j = 0; j < max_range; j++)
						for (k = 0; k < max_range; k++)
							M_1[i][j][k] = temperature_t(i,j,k,t,M_2);
				//отправляем направо
				if (rank < count - 1 && t < T_N)
				{
					MPI_Isend(&M_1[i_right-1][0][0], max_range*max_range, MPI_DOUBLE, rank+1, (t+1)*(rank*count+(rank+1)), MPI_COMM_WORLD, &req1);
				}
				//отправляем налево
				if (rank > 1 && t < T_N)
				{
					MPI_Isend(&M_1[i_left][0][0], max_range*max_range, MPI_DOUBLE, rank-1, (t+1)*(rank*count+(rank-1)), MPI_COMM_WORLD, &req2);
				}
				
			}
			if (ceil(t*delta_t) == 2)
				MPI_Send(&M_1[i_left][0][0],max_range*max_range*(i_right-i_left),MPI_DOUBLE,0,rank*2,MPI_COMM_WORLD);
			if (ceil(t*delta_t) == 4)
				MPI_Send(&M_1[i_left][0][0],max_range*max_range*(i_right-i_left),MPI_DOUBLE,0,rank*4,MPI_COMM_WORLD);
			if (ceil(t*delta_t) == 6)
				MPI_Send(&M_1[i_left][0][0],max_range*max_range*(i_right-i_left),MPI_DOUBLE,0,rank*6,MPI_COMM_WORLD);
			if (ceil(t*delta_t) == 10)
				MPI_Send(&M_1[i_left][0][0],max_range*max_range*(i_right-i_left),MPI_DOUBLE,0,rank*10,MPI_COMM_WORLD);

			//Пишем свой блок в M_2
			for (i = 0; i < max_range; i++)
				for (j = 0; j < max_range; j++)
					for (k = 0; k < max_range; k++)
						M_2[i][j][k] = M_1[i][j][k];
		}
		MPI_Send(&M_1[i_left][0][0],max_range*max_range*(i_right-i_left),MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
	}
	else if (rank == 0)
	{
		//print_matrix(M_2, t_end-t_start,file5);
		t_start = MPI_Wtime();
		for (i = 1; i < count; i++)
		{
		 	i_left = 0 + (max_range / work_count) * (i - 1);
		 	i_right = 0 + (max_range / work_count) * i;
			if (i == count - 1 && i_right < max_range)
				i_right = max_range;
		 	MPI_Recv(&M_2[i_left][0][0],max_range*max_range*(i_right-i_left),MPI_DOUBLE,i,i*2,MPI_COMM_WORLD,&status1);
		}
		print_temperature(M_2, 2, file1);
		for (i = 1; i < count; i++)
		{
		 	i_left = 0 + (max_range / work_count) * (i - 1);
		 	i_right = 0 + (max_range / work_count) * i;
			if (i == count - 1 && i_right < max_range)
				i_right = max_range;
		 	MPI_Recv(&M_2[i_left][0][0],max_range*max_range*(i_right-i_left),MPI_DOUBLE,i,i*4,MPI_COMM_WORLD,&status1);
		}
		print_temperature(M_2, 4, file2);
		for (i = 1; i < count; i++)
		{
		 	i_left = 0 + (max_range / work_count) * (i - 1);
		 	i_right = 0 + (max_range / work_count) * i;
			if (i == count - 1 && i_right < max_range)
				i_right = max_range;
		 	MPI_Recv(&M_2[i_left][0][0],max_range*max_range*(i_right-i_left),MPI_DOUBLE,i,i*6,MPI_COMM_WORLD,&status1);
		}
		print_temperature(M_2, 6, file3);
		for (i = 1; i < count; i++)
		{
		 	i_left = 0 + (max_range / work_count) * (i - 1);
		 	i_right = 0 + (max_range / work_count) * i;
			if (i == count - 1 && i_right < max_range)
				i_right = max_range;
		 	MPI_Recv(&M_2[i_left][0][0],max_range*max_range*(i_right-i_left),MPI_DOUBLE,i,i*10,MPI_COMM_WORLD,&status1);
		}
		print_temperature(M_2, 10, file4);
		for (i = 1; i < count; i++)
		{
		 	i_left = 0 + (max_range / work_count) * (i - 1);
		 	i_right = 0 + (max_range / work_count) * i;
			if (i == count - 1 && i_right < max_range)
				i_right = max_range;
		 	MPI_Recv(&M_2[i_left][0][0],max_range*max_range*(i_right-i_left),MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status1);
		}
		t_end = MPI_Wtime();
		print_matrix(M_2, t_end-t_start,file5);
		printf("Done \n");
		printf("spent time = %3f \n", t_end-t_start);
		
	}
	MPI_Finalize();
	return 0;
}

