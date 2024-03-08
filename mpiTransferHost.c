#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include "pm18MpiHeader.h"

#include "mpiTransferHost.h"

static MPI_Request     *requestSendMessages;
static MPI_Status       *statusSendMessages;
static MPI_Request     *requestRecvMessages;
static MPI_Status       *statusRecvMessages;

static MPI_Request *requestSendMessagesGrad;
static MPI_Status   *statusSendMessagesGrad;
static MPI_Request *requestRecvMessagesGrad;
static MPI_Status   *statusRecvMessagesGrad;

int InitHostsInterconnectRequestsMPI(MpiHostInterconnectType hostInterconnect, MPI_Comm comm)
{
	int i;

	// ÈÍÈÖÈÀËÈÇÀÖÈß ÑÒÐÓÊÒÓÐ ÇÀÏÓÑÊÀ È ÏÐÎÂÅÐÊÈ ÇÀÂÅÐØÅÍÍÎÑÒÈ ÇÀÏÐÎÑÎÂ ÍÀ ÏÅÐÅÄÀ×Ó ÄÀÍÍÛÕ
	if(hostInterconnect.sendMessagesNumber > 0)
	{
		requestSendMessages = (MPI_Request *)malloc(hostInterconnect.sendMessagesNumber * sizeof(MPI_Request)); if(requestSendMessages == NULL) exit(0);
		statusSendMessages  = (MPI_Status  *)malloc(hostInterconnect.sendMessagesNumber * sizeof(MPI_Status )); if(statusSendMessages  == NULL) exit(0);
		for(i=0; i<hostInterconnect.sendMessagesNumber; i++)
		{
			MPI_Send_init(hostInterconnect.sendBufferPointers[i], hostInterconnect.sendMessageLengths[i]*5, MPI_DOUBLE, 
				          hostInterconnect.sendMessageHostsIDs[i], 0, comm, 
						  &requestSendMessages[i]);		
		} // for i
	}
	else
	{
		requestSendMessages = NULL;
		statusSendMessages  = NULL;
	}

	// ÈÍÈÖÈÀËÈÇÀÖÈß ÑÒÐÓÊÒÓÐ ÇÀÏÓÑÊÀ È ÏÐÎÂÅÐÊÈ ÇÀÂÅÐØÅÍÍÎÑÒÈ ÇÀÏÐÎÑÎÂ ÍÀ ÏÐÈÅÌ ÄÀÍÍÛÕ
	if(hostInterconnect.recvMessagesNumber > 0)
	{
		requestRecvMessages = (MPI_Request *)malloc(hostInterconnect.recvMessagesNumber * sizeof(MPI_Request)); if(requestRecvMessages == NULL) exit(0);
		statusRecvMessages  = (MPI_Status  *)malloc(hostInterconnect.recvMessagesNumber * sizeof(MPI_Status )); if(statusRecvMessages  == NULL) exit(0);
		for(i=0; i<hostInterconnect.recvMessagesNumber; i++)
		{
			MPI_Recv_init(hostInterconnect.recvBufferPointers[i], hostInterconnect.recvMessageLengths[i]*5, MPI_DOUBLE, 
				          hostInterconnect.recvMessageHostsIDs[i], 0, comm, 
						  &requestRecvMessages[i]);		
		} // for i
	}
	else
	{
		requestRecvMessages = NULL;
		statusRecvMessages  = NULL;
	}

	return 0;
}

int HostsInterconnectStartWaitAll(MpiHostInterconnectType hostInterconnect, MPI_Comm comm)
{

	if(hostInterconnect.sendMessagesNumber > 0) MPI_Startall(hostInterconnect.sendMessagesNumber, requestSendMessages);
	if(hostInterconnect.recvMessagesNumber > 0) MPI_Startall(hostInterconnect.recvMessagesNumber, requestRecvMessages);
	if(hostInterconnect.sendMessagesNumber > 0) MPI_Waitall(hostInterconnect.sendMessagesNumber,  requestSendMessages, statusSendMessages);
	if(hostInterconnect.recvMessagesNumber > 0) MPI_Waitall(hostInterconnect.recvMessagesNumber,  requestRecvMessages, statusRecvMessages);

	return 0;
}

int RebootHostsInterconnectRecvRequestsMPI(MpiHostInterconnectType *hostInterconnect, double *QphRecv, MPI_Comm comm)
{
	int rankMPI;

	MPI_Comm_rank(comm, &rankMPI);

	// ÓÄÀËÅÍÈÅ ÈÄÅÍÒÈÔÈÊÀÒÎÐÎÂ ÑÒÀÐÛÕ ÇÀÏÐÎÑÎÂ
	if(hostInterconnect->recvMessagesNumber > 0)
	{
		int i;
		for(i=0; i<hostInterconnect->recvMessagesNumber; i++)
			MPI_Request_free(&requestRecvMessages[i]);
	}
	// ÓÄÀËÅÍÈÅ ÈÄÅÍÒÈÔÈÊÀÒÎÐÎÂ ÑÒÀÐÛÕ ÇÀÏÐÎÑÎÂ

	// ÈÇÌÅÍÅÍÈÅ ÀÄÐÅÑÎÂ Â ÁÓÔÅÐÀÕ ÏÐÈÅÌÀ ÄÀÍÍÛÕ
	if(hostInterconnect->recvMessagesNumber > 0)
	{
		int i, j;

		hostInterconnect->recvValuesBuffer = QphRecv;

		for(j=0, i=0; i<hostInterconnect->recvMessagesNumber; i++)
		{
			hostInterconnect->recvBufferPointers[i] = hostInterconnect->recvValuesBuffer + j * 5;
			j += hostInterconnect->recvMessageLengths[i];
		}

		if(j != hostInterconnect->recvBufferTotalLength)
		{ fprintf(stderr, "[%d]: ERROR hostInterconnect->recvBufferTotalLength\n", rankMPI); fflush(stderr); exit(0); }
	}
	// ÈÇÌÅÍÅÍÈÅ ÀÄÐÅÑÎÂ Â ÁÓÔÅÐÀÕ ÏÐÈÅÌÀ ÄÀÍÍÛÕ
	
	// ÈÍÈÖÈÀËÈÇÀÖÈß ÇÀÏÐÎÑÎÂ ÍÀ ÀÑÈÍÕÐÎÍÍÛÉ ÏÐÈÅÌ ÄÀÍÍÛÕ
	if(hostInterconnect->recvMessagesNumber > 0)
	{
		int i;

		for(i=0; i<hostInterconnect->recvMessagesNumber; i++)
		{
			MPI_Recv_init(hostInterconnect->recvBufferPointers[i], hostInterconnect->recvMessageLengths[i]*5, MPI_DOUBLE, 
				          hostInterconnect->recvMessageHostsIDs[i], 0, comm, 
						  &requestRecvMessages[i]);
		} // for i
	}
	// ÈÍÈÖÈÀËÈÇÀÖÈß ÇÀÏÐÎÑÎÂ ÍÀ ÀÑÈÍÕÐÎÍÍÛÉ ÏÐÈÅÌ ÄÀÍÍÛÕ

	MPI_Barrier(comm);
	if(rankMPI == 0) { printf("REBOOT RECV REQUESTS\n"); fflush(stdout); }
	MPI_Barrier(comm);	
	
	return 0;
}

int InitHostsInterconnectGradientRequestsMPI(MpiHostInterconnectType *hostInterconnect, double *GradRecv, MPI_Comm comm)
{
	int rankMPI;

	int i, j;

	MPI_Comm_rank(comm, &rankMPI);

	if(hostInterconnect->sendMessagesNumber > 0)
	{
		hostInterconnect->sendBufferPointersGrad = (void  **)malloc(hostInterconnect->sendMessagesNumber    *      sizeof(void *)); if(hostInterconnect->sendBufferPointersGrad == NULL) exit(0);
		hostInterconnect->sendValuesBufferGrad   = (double *)malloc(hostInterconnect->sendBufferTotalLength * 15 * sizeof(double)); if(hostInterconnect->sendValuesBufferGrad   == NULL) exit(0);
		for(i=0; i<hostInterconnect->sendBufferTotalLength * 15; i++) hostInterconnect->sendValuesBufferGrad[i] = 0.0; 

		for(j=0, i=0; i<hostInterconnect->sendMessagesNumber; i++)
		{
			hostInterconnect->sendBufferPointersGrad[i] = (void *)(&hostInterconnect->sendValuesBufferGrad[j * 15]);
			j += hostInterconnect->sendMessageLengths[i];
		} // for i

		requestSendMessagesGrad = (MPI_Request *)malloc(hostInterconnect->sendMessagesNumber * sizeof(MPI_Request)); if(requestSendMessagesGrad == NULL) exit(0);
		statusSendMessagesGrad  = (MPI_Status  *)malloc(hostInterconnect->sendMessagesNumber * sizeof(MPI_Status )); if(statusSendMessagesGrad  == NULL) exit(0);

		for(i=0; i<hostInterconnect->sendMessagesNumber; i++)
		{
			MPI_Send_init(hostInterconnect->sendBufferPointersGrad[i], hostInterconnect->sendMessageLengths[i]*15, MPI_DOUBLE, 
				          hostInterconnect->sendMessageHostsIDs[i], 257, comm, 
						  &requestSendMessagesGrad[i]);
		} // for i
	}
	else
	{
		hostInterconnect->sendBufferPointersGrad = NULL; 
		hostInterconnect->sendValuesBufferGrad = NULL;

		requestSendMessagesGrad = NULL;
		statusSendMessagesGrad  = NULL;
	}

	if(hostInterconnect->recvMessagesNumber > 0)
	{
		hostInterconnect->recvValuesBufferGrad = GradRecv;
		hostInterconnect->recvBufferPointersGrad = (void  **)malloc(hostInterconnect->recvMessagesNumber * sizeof(void *)); if(hostInterconnect->recvBufferPointersGrad == NULL) exit(0);

		for(j=0, i=0; i<hostInterconnect->recvMessagesNumber; i++)
		{
			hostInterconnect->recvBufferPointersGrad[i] = hostInterconnect->recvValuesBufferGrad + j * 15;
			j += hostInterconnect->recvMessageLengths[i];
		}		

		requestRecvMessagesGrad = (MPI_Request *)malloc(hostInterconnect->recvMessagesNumber * sizeof(MPI_Request)); if(requestRecvMessagesGrad == NULL) exit(0);
		statusRecvMessagesGrad  = (MPI_Status  *)malloc(hostInterconnect->recvMessagesNumber * sizeof(MPI_Status )); if(statusRecvMessagesGrad  == NULL) exit(0);

		for(i=0; i<hostInterconnect->recvMessagesNumber; i++)
		{
			MPI_Recv_init(hostInterconnect->recvBufferPointersGrad[i], hostInterconnect->recvMessageLengths[i]*15, MPI_DOUBLE, 
				          hostInterconnect->recvMessageHostsIDs[i], 257, comm, 
						  &requestRecvMessagesGrad[i]);
		} // for i
	}
	else
	{
		hostInterconnect->recvValuesBufferGrad = NULL;
		hostInterconnect->recvBufferPointersGrad = NULL;

		requestRecvMessagesGrad = NULL;
		statusRecvMessagesGrad  = NULL;
	}	

	MPI_Barrier(comm);
	if(rankMPI == 0) { printf("INIT GRADIENTS DATA TRANSFER REQUESTS\n"); fflush(stdout); }
	MPI_Barrier(comm);	
	
	return 0;
}

int HostsInterconnectGradientsStartWaitAll(MpiHostInterconnectType hostInterconnect, MPI_Comm comm)
{

	if(hostInterconnect.sendMessagesNumber > 0) MPI_Startall(hostInterconnect.sendMessagesNumber, requestSendMessagesGrad);
	if(hostInterconnect.recvMessagesNumber > 0) MPI_Startall(hostInterconnect.recvMessagesNumber, requestRecvMessagesGrad);
	if(hostInterconnect.sendMessagesNumber > 0) MPI_Waitall(hostInterconnect.sendMessagesNumber,  requestSendMessagesGrad, statusSendMessagesGrad);
	if(hostInterconnect.recvMessagesNumber > 0) MPI_Waitall(hostInterconnect.recvMessagesNumber,  requestRecvMessagesGrad, statusRecvMessagesGrad);

	return 0;
}
