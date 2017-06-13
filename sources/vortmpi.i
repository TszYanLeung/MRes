
      integer nproc       ! number of processors
      integer myid       ! node id
      integer root       ! root node id
      integer comm       ! mpi global communicator 
      integer real_t     ! real mpi data type
      integer realsp_t   ! real single precision formovies mpi data type
      integer int_t      ! integer mpi data type
      integer char_t     ! character mpi data type
      integer pack_t     ! packed mpi data type
      integer bool_t     ! logical mpi data type
      integer stat_s     ! status size
      integer sum_oper     ! sum reduce operation
      integer max_oper     ! max reduce operation
      common /ptool0/ nproc,myid,root,comm
      common /ptool1/ real_t,realsp_t,int_t,char_t,pack_t,bool_t,stat_s
      common /ptool2/ sum_oper,max_oper
