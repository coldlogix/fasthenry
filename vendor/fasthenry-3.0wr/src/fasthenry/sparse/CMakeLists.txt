cmake_minimum_required(VERSION 2.6)

project(sparse C)


set (HFILES 
	spConfig.h spDefs.h spMatrix.h
)

set (CFILES
	spAllocate.c spBuild.c spFactor.c spOutput.c spSolve.c spUtils.c
	 spCompat.c spFortran.c
)

add_library(sparse
	${CFILES} ${HFILES}
)

