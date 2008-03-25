/* $Header: /home/cvs/bp/oofem/oofemlib/src/femcmpnn.C,v 1.12.4.1 2004/04/05 15:19:43 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

/*
 The original idea for this class comes from 
  Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 PhD Thesis, EPFL, Lausanne, 1992.
*/

//   file FEMCMPNN.C

#include "femcmpnn.h"
#include "domain.h"
#include "cltypes.h"
#include "engngm.h"
#include "logger.h"
#include "datastream.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#endif



contextIOResultType
FEMComponent::saveContext (DataStream* stream, ContextMode mode, void *obj)
{ 
 if (stream == NULL) 
  THROW_CIOERR(CIO_IOERR);
 int type_id = this->giveClassID ();
 // write class header
 if (!stream->write(&type_id,1)) THROW_CIOERR(CIO_IOERR);

 if (mode & CM_Definition ) {
   if (!stream->write(&number,1)) THROW_CIOERR(CIO_IOERR);
 }
 return CIO_OK;
}

contextIOResultType   
FEMComponent::restoreContext(DataStream* stream, ContextMode mode, void *obj)
{
 int class_id;
 int type_id = this->giveClassID ();
 // read class header
 if (!stream->read(&class_id,1)) THROW_CIOERR(CIO_IOERR);
 if (class_id != type_id) 
  THROW_CIOERR(CIO_BADVERSION);

 if (mode & CM_Definition ) {
   if (!stream->read(&number,1)) THROW_CIOERR(CIO_IOERR);
 }
 
 
 return CIO_OK;
}


int
FEMComponent :: giveInputRecordString(std::string &str, bool keyword)
{
 char buffer[512];

 if(keyword == true){
  sprintf(buffer, "%s %d", this -> giveInputRecordName(), this -> giveNumber());
  str = buffer;
 }

 return 1;
}



void 
FEMComponent :: error (const char* file, int line, const char *format, ...) const
{
  char buffer[MAX_ERROR_MSG_LENGTH];
	va_list args;

	va_start(args, format);
	vsprintf(buffer, format, args);
	va_end(args);
  
  __OOFEM_ERROR4 (file, line, "Class: %s, number: %d\n%s",giveClassName(),giveNumber(), buffer);
}

void 
FEMComponent :: warning (const char* file, int line, const char *format, ...) const
{
  char buffer[MAX_ERROR_MSG_LENGTH];
	va_list args;

	va_start(args, format);
	vsprintf(buffer, format, args);
	va_end(args);
  
  __OOFEM_WARNING4 (file, line, "Class: %s, number: %d\n%s",giveClassName(),giveNumber(), buffer);
}

/*
// ___DO_NOT_CHANGE_FOLLOWING_LINES
// they contain typeInfo informations for base class
const TypeInfo  FEMComponent::infoFEMComponent ("FEMComponent",0,NULL);

const TypeInfo* FEMComponent::getTypeInfo() const 
 {return &FEMComponent::infoFEMComponent;}
// end of typeInfo section
*/
