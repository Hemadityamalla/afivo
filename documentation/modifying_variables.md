# Modyfing variables

# Introduction

Afivo supports cell-centered and face-centered variables, which can be added with `m_af_core::af_add_cc_variable()` and `m_af_core::af_add_fc_variable()`.

# Changing variables on a single box

A single box contains an array `cc` with cell-centered variables and an array
`fc` with face centered variables, see the definition in `m_af_types::box_t`.
These variables can for example be modified as shown below:

```{f90}
    subroutine change_values(box)
       type(box2_t), intent(inout) :: box
       integer                     :: i, j

       ! Loop over the cells
       do j = 1, box%n_cell
          do i = 1, box%n_cell
             ! Set cell-centered variable 1
             box%cc(i, j, 1) = i + j + 1.0

             ! Set cell-centered variable 2
             box%cc(i, j, 2) = i - j - 1.0

             ! Set face-centered variable 1 in the x-direction
             box%fc(i, j, 1, 1) = 1.0

             ! Set face-centered variable 1 in the y-direction
             box%fc(i, j, 2, 1) = 0.0
          end do
       end do
    end subroutine change_values
```

# Calling a routine on all boxes of the tree

If a user has written a routine as shown above, then this routine can be passed
to `m_af_utils::af_loop_box`:

    call af_loop_box(tree, change_values)

The routine will then be called on all the boxes of the tree. See `m_af_utils`
for alternatives to `af_loop_box`.

# Manually looping over all the boxes in a tree

A user can also manually loop over all the boxes in a tree to perform some
operations. A user can for example loop over all the leaves (boxes without
children) like this:

```{f90}
    do lvl = 1, tree%highest_lvl
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call change_values(tree%boxes(id))
       end do
    end do
```

For simplicity, OpenMP commands
for [parallelization](@ref documentation/parallelization.md) have been omitted
in this example.
