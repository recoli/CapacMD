/*****************************************************************************
 *  This file is part of the CapacMD program.                                *  
 *  Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>                         *  
 *                                                                           *
 *  Filename:  my_malloc.c                                                   *
 *  Function:  allocate memory with error checking                           *
 *  Version:   1.0                                                           *
 *  Updated:   2015-Jun-30                                                   *
 *  License:   GNU Public License, version 2                                 *
 *                                                                           *  
 *  This program is free software; you can redistribute it and/or modify     *  
 *  it under the terms of the GNU General Public License as published by     *  
 *  the Free Software Foundation; either version 2 of the License, or        *  
 *  (at your option) any later version.                                      *  
 *                                                                           *
 *  This program is distributed in the hope that it will be useful,          *  
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *  
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *  
 *  GNU General Public License for more details.                             *  
 *                                                                           *
 *  You should have received a copy of the GNU General Public License along  *  
 *  with this program; if not, write to the Free Software Foundation, Inc.,  *  
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.              *  
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

//=============================
// my own malloc function
//=============================

void* my_malloc(size_t bytes) 
{
	void* ptr = malloc(bytes);
	if(ptr == NULL) 
	{
		printf ("Error: could not allocate memory!\n");
		exit(1);
	} 
	else 
	{
		return ptr;
	}
}

void* my_malloc_2(size_t bytes, char *word)
{
	void* ptr = malloc(bytes);

#ifdef DEBUG
	printf("size of alloc (%s) = %zu MB\n", word, bytes / 1000000);
#endif

	if(ptr == NULL) 
	{
		printf ("Error: could not allocate memory for %s !\n", word);
		exit(1);
	} 
	else 
	{
		return ptr;
	}
}

