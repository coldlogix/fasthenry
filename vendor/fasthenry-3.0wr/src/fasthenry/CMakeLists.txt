cmake_minimum_required(VERSION 2.6)

project(fasthenry C)

# These CFLAGS should be fine for anything other than MinGW.
#CFLAGS = -O -g -DFOUR
add_definitions(-DFOUR)

add_subdirectory(sparse)

# For MinGW, there is no getrusage, and we'll skip a gethostname call
# to avoid having to link with winsock.
#CFLAGS = -O -g -DFOUR -DNO_RUSAGE -DNO_GETHOSTNAME

set (HEADER 
	induct.h cmplx.h resusage.h
)

set (MULHEAD 
	mulStruct.h mulGlobal.h patran.h resusage.h
)

set (NONUNIHEAD 
	gp.h
)

set (SRC
	induct.c gmres.c savemat_mod.c readGeom.c joelself.c writefastcap.c 
       SetupMulti.c calcp.c SetupComputePsi.c mulSetup.c BreakupSeg.c 
       Precond.c addgroundplane.c findpaths.c fillM.c cx_ludecomp.c 
       parse_command_line.c default_opts.c Prec_cost.c hole.c 
       dist_betw_fils.c mutual.c newPrecond.c deg_mutual.c 
       barnoldi.c regurgitate.c
)

set (MSRC
	mulGlobal.c mulDo.c mulMulti.c mulLocal.c mulMats.c direct.c
        uglieralloc.c capsolve.c
)

set (NONUNISRC
	find_nonuni_path.c read_tree.c contact.c
)

add_executable(fasthenry
    ${SRC} ${MSRC} ${NONUNISRC}
)

target_link_libraries(fasthenry sparse)


