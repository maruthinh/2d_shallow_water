#declare the variable
CC=g++

CFLAGS=-c -Wall

all: shallow2d

shallow2d: initialize_variables.o read_grid.o grid_computations.o circular_dam_break.o ic_oblique_hydraulic_jump.o ic_dam_break.o bc_prescribed_inflow.o dependent_variables.o bc_transmitive.o bc_wall.o bc_circular_dam_break.o bc_oblique_hydraulic_jump.o bc_dam_break.o avg_conv_flux1.o time_step.o var_diffs.o avg_conv_flux1.o avg_conv_flux2.o diss_movers.o diss_movers_h1.o diss_movers_le1.o movers_h1.o diss_movers_h2_prim.o diss_movers_h2_conv.o diss_movers_le2_conv.o diss_movers_le2_prim.o flux_boundary.o solver_explicit.o solver.o write_solution.o main.o 
	$(CC) initialize_variables.o read_grid.o grid_computations.o circular_dam_break.o ic_oblique_hydraulic_jump.o ic_dam_break.o bc_prescribed_inflow.o dependent_variables.o bc_transmitive.o bc_wall.o bc_circular_dam_break.o bc_oblique_hydraulic_jump.o bc_dam_break.o time_step.o var_diffs.o avg_conv_flux1.o avg_conv_flux2.o diss_movers_h1.o diss_movers.o diss_movers_le1.o movers_h1.o diss_movers_h2_prim.o diss_movers_h2_conv.o diss_movers_le2_conv.o diss_movers_le2_prim.o flux_boundary.o solver_explicit.o solver.o write_solution.o main.o -o shallow2d

initialize_variables.o:initialize_variables.cpp
	$(CC) $(CFLAGS) initialize_variables.cpp

read_grid.o:read_grid.cpp
	$(CC) $(CFLAGS) read_grid.cpp
	
grid_computations.o:grid_computations.cpp
	$(CC) $(CFLAGS) grid_computations.cpp

circular_dam_break.o:circular_dam_break.cpp
	$(CC) $(CFLAGS) circular_dam_break.cpp	

ic_oblique_hydraulic_jump.o:ic_oblique_hydraulic_jump.cpp
	$(CC) $(CFLAGS) ic_oblique_hydraulic_jump.cpp

ic_dam_break.o:ic_dam_break.cpp
	$(CC) $(CFLAGS) ic_dam_break.cpp
	
bc_prescribed_inflow.o:bc_prescribed_inflow.cpp
	$(CC) $(CFLAGS) bc_prescribed_inflow.cpp
	
dependent_variables.o:dependent_variables.cpp
	$(CC) $(CFLAGS) dependent_variables.cpp

bc_transmitive.o:bc_transmitive.cpp
	$(CC) $(CFLAGS) bc_transmitive.cpp
	
bc_wall.o:bc_wall.cpp
	$(CC) $(CFLAGS) bc_wall.cpp
		
bc_circular_dam_break.o:bc_circular_dam_break.cpp
	$(CC) $(CFLAGS) bc_circular_dam_break.cpp
	
bc_oblique_hydraulic_jump.o:bc_oblique_hydraulic_jump.cpp
	$(CC) $(CFLAGS) bc_oblique_hydraulic_jump.cpp
	
bc_dam_break.o:bc_dam_break.cpp
	$(CC) $(CFLAGS) bc_dam_break.cpp
	
time_step.o:time_step.cpp
	$(CC) $(CFLAGS) time_step.cpp
	
var_diffs.o:var_diffs.cpp
	$(CC) $(CFLAGS) var_diffs.cpp
	
avg_conv_flux1.o:avg_conv_flux1.cpp
	$(CC) $(CFLAGS) avg_conv_flux1.cpp
	
avg_conv_flux2.o:avg_conv_flux2.cpp
	$(CC) $(CFLAGS) avg_conv_flux2.cpp

diss_movers.o:diss_movers.cpp
	$(CC) $(CFLAGS) diss_movers.cpp
	
diss_movers_h1.o:diss_movers_h1.cpp
	$(CC) $(CFLAGS) diss_movers_h1.cpp	
	
movers_h1.o:movers_h1.cpp
	$(CC) $(CFLAGS) movers_h1.cpp	
	
diss_movers_le1.o:diss_movers_le1.cpp
	$(CC) $(CFLAGS) diss_movers_le1.cpp

diss_movers_h2_prim.o:diss_movers_h2_prim.cpp
	$(CC) $(CFLAGS) diss_movers_h2_prim.cpp	
	
diss_movers_h2_conv.o:diss_movers_h2_conv.cpp
	$(CC) $(CFLAGS) diss_movers_h2_conv.cpp			

diss_movers_le2_prim.o:diss_movers_le2_prim.cpp
	$(CC) $(CFLAGS) diss_movers_le2_prim.cpp	
	
diss_movers_le2_conv.o:diss_movers_le2_conv.cpp
	$(CC) $(CFLAGS) diss_movers_le2_conv.cpp			
	
flux_boundary.o:flux_boundary.cpp
	$(CC) $(CFLAGS) flux_boundary.cpp

solver_explicit.o:solver_explicit.cpp
	$(CC) $(CFLAGS) solver_explicit.cpp
	
solver.o:solver.cpp
	$(CC) $(CFLAGS) solver.cpp
	
write_solution.o:write_solution.cpp
	$(CC) $(CFLAGS) write_solution.cpp
	
main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

clean: 
	rm -rf *o shallow2d
