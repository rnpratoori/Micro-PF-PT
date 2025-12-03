// save_checkpoint
template <int dim> void Solid<dim>::save_checkpoint() {
  pcout << "Saving checkpoint..." << std::endl;

  // 1. Save triangulation and DoF handlers
  parallel::distributed::SolutionTransfer<dim, vectorType> sol_trans(
      dof_handler);
  parallel::distributed::SolutionTransfer<dim, vectorType> sol_trans_c(
      dof_handler_c);

  triangulation.prepare_coarsening_and_refinement();
  sol_trans.prepare_for_coarsening_and_refinement(solution_u);
  sol_trans_c.prepare_for_coarsening_and_refinement(
      solution_c1); // Saving c1, others similar?
  // Note: Full version saves c1, c2, c3. We need to handle all.

  // Actually, for simple checkpointing without refinement, we just save the
  // vectors. But deal.II checkpointing usually involves saving the
  // triangulation.

  // Let's follow the full version implementation if available.
  // Full version uses:
  /*
    triangulation.save("restart/mesh.triangulation");
    dof_handler.save("restart/dof_handler");
    // ... serialization of vectors ...
  */

  // Since we don't have the full version implementation handy in the context
  // (it was in microPF.cc which I viewed partially), I will implement a
  // standard deal.II checkpointing using serialization.

  std::string filename = output_directory + "/restart";
  std::ofstream out(filename.c_str());
  boost::archive::binary_oarchive oa(out);

  oa & triangulation;
  oa & dof_handler;
  oa & dof_handler_c;

  // Save time
  double current_time = time.current();
  double time_step = time.get_timestep();
  oa & current_time;
  oa & time_step;

  // Save solutions
  solution_u.print(std::cout); // Debug
  // Trilinos vectors might not be directly serializable with boost.
  // We usually save them using their own methods or copy to std::vector.
  // For distributed vectors, it's trickier.

  // Let's check how full version does it.
  // I need to view full/microPF.cc again to see save_checkpoint implementation.
}
