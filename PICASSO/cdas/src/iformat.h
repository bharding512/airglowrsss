//=============================================================================
// Module name : iformat.h
// Copyright (C) 1999 - 2003 Jan Soldan
//               251 65 Ondrejov - 236
//               Czech Republic
//=============================================================================
#ifndef _IFORMAT_H_
#define _IFORMAT_H_

#include "sbigdef.h"
#include "status.h"
//=============================================================================
class IFormat : public Status
{
 public:
 IFormat();
 virtual ~IFormat();

 virtual int Load(const char *pathName, int extension = 0) = 0;
 virtual int Save(const char *pathName) = 0;
};
//=============================================================================
#endif //_IFORMAT_H_
