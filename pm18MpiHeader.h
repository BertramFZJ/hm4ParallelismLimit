#ifndef _MpiHostInterconnectType_
typedef struct
{
	int hostsNumber;
	int hostID;	

	int recvBufferTotalLength, *RecvX, *RecvAGlobalIndex;
	int sendBufferTotalLength, *SendX, *SendAGlobalIndex, *SendALocalIndex;

	int    sendMessagesNumber;
    int  *sendMessageHostsIDs;
	int   *sendMessageLengths;
	void **sendBufferPointers;
	double  *sendValuesBuffer;
	
	int    recvMessagesNumber;
	int  *recvMessageHostsIDs;
	int   *recvMessageLengths;
	void **recvBufferPointers;
	double  *recvValuesBuffer;

	void **sendBufferPointersGrad;
	double  *sendValuesBufferGrad;
	void **recvBufferPointersGrad;
	double *recvValuesBufferGrad;
} MpiHostInterconnectType;
#define _MpiHostInterconnectType_
#endif

#ifndef _pm18MpiTopologyType_
typedef struct
{
	int rankWORLD;
	int sizeWORLD;
	
	int rankTASK;
	int sizeTASK;
	MPI_Comm commTASK;

	MpiHostInterconnectType hostInterconnect;
} pm18MpiTopologyType;
#define _pm18MpiTopologyType_
#endif
