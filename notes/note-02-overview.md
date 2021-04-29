A small project to do a pure port of ARPACK Symmetric solver to Julia - Part 2 - Goals

See [Overview](https://dgleich.micro.blog/2021/04/21/a-small-project.html) for a
high-level picture.

Here is a quick call graph of the various ARPACK functions for the symmetric
eigenvalue/eigenvector problem `dsaupd`

- `d` = double precision
- `s` = symmetric
- `aupd` = Arnoldi Update (since we'll call this in a loop to update the
internal Arnoldi information.)

> Note to the pedantic. Yes, I know the Arnoldi of a symmetric matrix is usually
> called the Lanczos process. So I should be saying Lanczos everywhere here. Except,
> maybe I really shouldn't be saying either and we could use the phrase:
> iteratively orthogonalized power basis instead of Arnoldi or
> three-term recurrent power basis instead of Lanczos? Either way, since the
> package is called ARPACK and I'm writing about that, I'll stick with
> Arnoldi here.

Back to the call graph.

    dsaupd calls
      dsaup2 calls
        dgetv0 calls printing, reverse communication
        dsaitr calls printing, blas, lapack
        dsapps calls printing, blas, lapack
        dsconv
        dseigt
        dsgets
        dsortr

        ivout - write out an integer vector
        arscnd - internal wrapper avoiding LAPACK second function
        dvout - write out a double vector

then we also need to handle getting the eigenvector information out of the computation.

    dseupd calls
      dsesrt
      dsortr - same as above

This suggests a plan.

Port the low-level/internal routines first. Then port the high-level
routines.

Wait what you say? Why not do it the opposite way? Well, I started doing it
that way and realized it's just harder ;). So we are going to switch
and do the easy stuff first!

Just as a quick overview of what the code actually does.

    dsaupd - high level interface
      dsaup2 - medium level interface, manages the overall process
        dgetv0 - gets the starting vector, or returns control to user to get it.
        dsaitr - haven't looked yet...
        dsapps - haven't looked yet...
        dsconv - check what has converged (easiest routine!)
        dseigt - haven't looked yet...
        dsgets - get shifts based on eigenvalue estimates
        dsortr - sorts a vector and/or applies the sort to an aux vector too

Also, let me note that the
[ARPACK documentation is stellar. Stellar.](https://twitter.com/dgleich/status/1382369397612765186)
So a few clicks and I could easily figure out what those other functions are doing.
