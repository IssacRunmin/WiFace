/*
 * (c) 2008-2011 Daniel Halperin <dhalperi@cs.washington.edu>
 */
#include "iwl_connector.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <signal.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <linux/netlink.h>
#include <sys/time.h>
#include <fcntl.h>
#define MAX_PAYLOAD 2048
#define SLOW_MSG_CNT 1
#define COMMUNICATE_FILE "/home/jennygroup/Matlab/com"
#define PING_COMMAND "nice -n 10 sudo ping 192.168.1.1 -r -i 0"

int sock_fd = -1;							// the socket
FILE* file = NULL;
FILE* communicate = NULL;
FILE* in = NULL;
int count = 0;
pthread_t id1;
char file_name[100];

void check_usage(int argc, char** argv);

FILE* open_file(char* filename, char* spec);
void *thread_ping();
void caught_signal(int sig);

void exit_program(int code);
void exit_program_err(int code, char* func);

int main(int argc, char** argv)
{
	/* Local variables */
	struct sockaddr_nl proc_addr, kern_addr;	// addrs for recv, send, bind
	struct cn_msg *cmsg;
	char buf[40960];
	char timebuf[128];
	int out;
	int ret,ret_ping;
	unsigned short l, l2;
	unsigned short int flag = 1;
	long int second, usecond;
	struct timeval start, buffer_time, write_time, interval;
	
	
	count = 0;

	/* Make sure usage is correct */
	check_usage(argc, argv);

	/* Open and write log file */
	strcpy(file_name,argv[1]);
	file = open_file(argv[1], "w");
	fclose(file);
	communicate = open_file(COMMUNICATE_FILE,"w");
	fprintf(communicate,"%s 0",argv[1]);
	fclose(communicate);
	out = open(argv[1], O_WRONLY|O_CREAT,S_IRUSR);
	//in = open_file("workspace/matlab/com","r");
	/* Setup the socket */
	sock_fd = socket(PF_NETLINK, SOCK_DGRAM, NETLINK_CONNECTOR);
	if (sock_fd == -1)
		exit_program_err(-1, "socket");

	/* Initialize the address structs */
	memset(&proc_addr, 0, sizeof(struct sockaddr_nl));
	proc_addr.nl_family = AF_NETLINK;
	proc_addr.nl_pid = getpid();			// this process' PID
	proc_addr.nl_groups = CN_IDX_IWLAGN;
	memset(&kern_addr, 0, sizeof(struct sockaddr_nl));
	kern_addr.nl_family = AF_NETLINK;
	kern_addr.nl_pid = 0;					// kernel
	kern_addr.nl_groups = CN_IDX_IWLAGN;
	
	/* Now bind the socket */
	if (bind(sock_fd, (struct sockaddr *)&proc_addr, sizeof(struct sockaddr_nl)) == -1)
		exit_program_err(-1, "bind");

	/* And subscribe to netlink group */
	{
		int on = proc_addr.nl_groups;
		ret = setsockopt(sock_fd, 270, NETLINK_ADD_MEMBERSHIP, &on, sizeof(on));
		if (ret)
			exit_program_err(-1, "setsockopt");
	}

	/* Set up the "caught_signal" function as this program's sig handler */
	signal(SIGINT, caught_signal);

	/* Poll socket forever waiting for a message */
	printf("Preprocess Done!");
	/*printf("Press Enter to start Collecting Data...\n");
	getchar();
	ret_ping = pthread_create(&id1, NULL, (void *)thread_ping,NULL);  
    if(ret_ping)  
    {  
        printf("Create pthread error!Please Ping manully...\n");  
    }  */
	communicate = open_file(COMMUNICATE_FILE,"w");
	fprintf(communicate,"%s 1",argv[1]);
	fclose(communicate);
	while (1)
	{
		//gettimeofday(&start, NULL);
		/* Receive from socket with infinite timeout */
		ret = recv(sock_fd, buf, sizeof(buf), 0);
		if (ret == -1){
			if (count == 5)
				exit_program_err(-1, "recv");
			else{
				count++;
				fprintf(stderr, "nobuffer %d times\n", count);
			}
		}
		else
			count = 0;
		/* Pull out the message portion and print some stats */
		cmsg = NLMSG_DATA(buf);
		//if (count % SLOW_MSG_CNT == 0)
			//printf("received %d bytes: id: %d val: %d seq: %d clen: %d\n", cmsg->len, cmsg->id.idx, cmsg->id.val, cmsg->seq, cmsg->len);
		/* Log the data to file */
		l = (unsigned short) cmsg->len;
		l2 = htons(l);
		/*fwrite(&l2, 1, sizeof(unsigned short), out);
		ret = fwrite(cmsg->data, 1, l, out);*/
		write(out, &l2, sizeof(l2));
		ret = write(out, &cmsg->data, l);
		//if (count % 100 == 0)
			//printf("wrote %d bytes [msgcnt=%u]\n", ret, count);
		//++count;
		if (ret != l)
			exit_program_err(1, "fwrite");
	}
	exit_program(0);
	return 0;
}

void *thread_ping(void)
{
	system(PING_COMMAND);
}

void check_usage(int argc, char** argv)
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s <output_file>\n", argv[0]);
		exit_program(1);
	}
}

FILE* open_file(char* filename, char* spec)
{
	FILE* fp = fopen(filename, spec);
	if (!fp)
	{
		perror("fopen");
		exit_program(1);
	}
	return fp;
}

void caught_signal(int sig)
{
	//fprintf(stderr, "msgcnt=%u\n", count);
	fprintf(stderr, "Caught signal %d\n", sig);
	//pthread_kill(&id1,0);
	exit_program(0);
}

void exit_program(int code)
{

	if(in)
	{
		fclose(in);
		in = NULL;
	}
	if (sock_fd != -1)
	{
		close(sock_fd);
		sock_fd = -1;
	}
	communicate = open_file(COMMUNICATE_FILE,"w");
	fprintf(communicate,"%s 2",file_name);
	fclose(communicate);
	
	exit(code);
}

void exit_program_err(int code, char* func)
{
	perror(func);
	exit_program(code);
}
