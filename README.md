# loads
Code from *Effects of dynamic power electronic load models for power systems analysis: A case study with ZIP-E loads*.
## Interactive plots:
- [Figure 2](https://reid.xz.ax/loads/transient2.html): Transient simulation for $dynpi$ lines at $load\\:scale\in\\{0.2, 0.5, 0.8\\}$.
- [Figure 3](https://reid.xz.ax/loads/transient1.html): Transient simulation for both $statpi$ and $dynpi$ lines. Includes slider for $load\\:scale$.
- [Figure 4](https://reid.xz.ax/loads/eigenvalues2.html): Eigenvalues for $dynpi$ lines at $load\\:scale\in\\{0.2, 0.5, 0.8\\}$.
- [Figure 5](https://reid.xz.ax/loads/eigenvalues1.html): Eigenvalues for both $statpi$ and $dynpi$ lines. Includes slider for $load\\:scale$.

## Reproducing our results

There are two steps: running the simulations and plotting the results.

First, run the simulations. This may take a while or a long while depending on how many cores your computer has - it took about ADD TIME HERE on our 128-core server. `-t $(nproc)` just ensures that julia uses as many threads as you have cores. 
```bash
julia -t $(nproc) experiments/9bus_sims.jl
```

Once you've run the simulations, results will be stored in the data folder. You can now read them and make plots.
```bash
julia experiments/9bus_results.jl
```

This will save the plots as `.html` files which you can open in your browser.

You'll notice that there's some additional data in these plots - there is a second column or row of plots containing the same results, but for the case with all synchronous machines. If you don't want this, you can set `just_like_paper=true` on line 20 of `9bus_results.jl`.
