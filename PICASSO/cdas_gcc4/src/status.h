//=============================================================================
// Module name : status.h
// Copyright (C) 2003 Jan Soldan
//               251 65 Ondrejov - 236
//               Czech Republic
//=============================================================================
#ifndef _STATUS_H_
#define _STATUS_H_

#include <base.h>
//=============================================================================
typedef enum{
									 
	S_NO_ERROR,
	S_PARAM_ERR,
	S_OPER_ERR,

	S_OBJ_ALLOC_ERR,
	S_OBJ_SIZE_ERR,
	S_OBJ_COPY_ERR,
	S_OBJ_DEF_ERR,
	S_OBJ_TRUNC_ERR,

	S_FILE_OPEN_ERR,
	S_FILE_READ_ERR,
	S_FILE_FORMAT_ERR,
	S_MAX_ITER_ERR,

	S_KEYWORD_ERR,
	S_VALUE_ERR,
	S_UNSUPPORTED_ERR,
	S_CANCEL

}STATUS_CODE;
//=============================================================================
// Status class
//=============================================================================
class Status
{
 protected:
 int                    m_status;

 public:
	                Status(){m_status = S_NO_ERROR;}
 virtual               ~Status(){}

 inline int             GetStatus(){return(m_status);}
 inline void            SetStatus(int status){m_status = status;}
 inline char           *GetStatusString();
};
//=============================================================================
#endif // _STATUS_H_
