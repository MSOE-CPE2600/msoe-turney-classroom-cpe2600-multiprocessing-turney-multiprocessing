/******************************************************************************
 * File: multi_mandel.c
 * Author: Zane Rothe
 * Description: multiprocessing of mandelbrot image movie (50 frames)
 * Arguments: enter number of desired concurrent processes
 * Class: CPE 2600-111
 * Assignment: Lab 11
 * Date: 11/20/24
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <semaphore.h>
#include <time.h>
#include <math.h>


int main(int argc, char* argv[])
{
	sem_t *sem;
    sem_unlink("/mysem"); //unlink any existing semaphore
    sem = sem_open("/mysem", O_CREAT, 0666, atoi(argv[1])); //create semaphore 
    if (sem == SEM_FAILED)
    {
        perror("sem_open");
        return 1;
    }

	int* frames= mmap(0, sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
	*frames=0; // stores current number of frames

	int orig_pid=getpid(); //remember initial parent process
	while (*frames<(50))
	{
		sem_wait(sem); // wait for available processer
		int pid=fork(); //make new process
		if(pid==0) // child only
		{
			(*frames)++; 
			float scale=6.7*exp(-0.18*(double)(*frames)); // change scale
			char scale_str[8];
			sprintf(scale_str, "%.5f", scale);
			char file_str[16];
			sprintf(file_str, "mandel%d.jpg", *frames); // change outfile
			//printf("processing %d\n",getpid());
			usleep(100000);
			execlp("./mandel","./mandel","-x","-0.5622","-y","0.6429","-s",scale_str,"-o",file_str,NULL);
			//execlp("./test","test",file_str,NULL); //test process
		}
		else
		{
			usleep(100000);
			//printf("waiting %d\n",getpid());
		}
	}
	if (getpid()!=orig_pid) // end any outstanding processes besides parent
	{
		kill(getpid(),SIGKILL);
	}
	
	int semval;
	sem_getvalue(sem,&semval);
	while(semval<atoi(argv[1])) //wait until semaphores back to initial value
	{	
		//printf("%d\n",semval);
		sleep(1); //do nothing
		sem_getvalue(sem,&semval);
	}
	printf("\ndone\n");

	munmap(frames, sizeof(int));
    sem_close(sem);
    sem_unlink("/mysem");
    return 0;
}