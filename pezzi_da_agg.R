



Again, similar structure:

```{r, fig.height=7, fig.width=7}
bimodal_emulators <- bimodal_emulator_from_data(all_output, output_names, ranges)
```

The function determines the bimodal structure across the output space, working
out the proportion of points in each mode. Then it looks at each output on its own,
and figures out if the output is bimodal. For the ones it doesn't think are bimodal
it trains stochastic emulators to them; for the ones it does think are bimodal it
trains an emulator to each mode. It includes unimodal emulators in both the modes,
so that we can treat a mode in isolation.

We end up with three objects: mode1, mode2, and prop. mode1 and mode2 each contain a
set of mean and variance emulators. Prop is a single emulator (determininstic) that
emulates the proportion in each mode across the space.

Quick look at the actives: note that the later time points have different active variables!
  
  ```{r, fig.height=7, fig.width=6}
plot_actives(bimodal_emulators$mode1$expectation)
plot_actives(bimodal_emulators$mode2$expectation)
```

```{r, fig.height=7, fig.width=9}
plot(bimodal_emulators$mode1$expectation$R350, params = c('alpha', 'sigma'))
plot(bimodal_emulators$mode2$expectation$R350, params = c('alpha', 'sigma'))
```

Bigggggg difference between the emulators.

Let's have a quick look at the prop emulator: it is just a bog-standard emulator.

```{r, fig.height=7, fig.width=9}
emulator_plot(bimodal_emulators$prop, params = c('alpha', 'sigma'))
```

At points where alpha (disease death) and gamma (recovery) are high, the first mode
dominates: this is the mode where the disease dies out. But there is always a
non-negligible risk of the disease dying out, since the proportion never drops
below 0.4

Might as well do validation diagnostics to point something out.
Also using subset_emulators to remove the I450 and R450 outputs.

```{r}
validation_diagnostics(subset_emulators(bimodal_emulators, output_names[-c(7,14)]), targets, all_valid)
```

You may note that we seem to have more points plotted than we have validation points.
This is because we do diagnostics on EACH mode for each output. In fact, for I350 (in
my running of it) you can actually see the mode separation in the plots and frequently you
can see the different emulators by the step change in uncertainty.

Let's finish off with some plots, akin to what we did earlier:
  
  ```{r, fig.height=7, fig.width=7}
em_preds1 <- purrr::map_dbl(bimodal_emulators$mode1$expectation, ~.$get_exp(chosen_df))
em_vars1 <- purrr::map_dbl(bimodal_emulators$mode1$expectation, ~.$get_cov(chosen_df))
em_preds2 <- purrr::map_dbl(bimodal_emulators$mode2$expectation, ~.$get_exp(chosen_df))
em_vars2 <- purrr::map_dbl(bimodal_emulators$mode2$expectation, ~.$get_cov(chosen_df))
plot(0:500, ylim=c(0,700), ty="n", xlab = "Time", ylab = "Number")
for(j in 3:4) for(i in 1:100) lines(0:500, solution[,j,i], col=(3:4)[j-2], lwd=0.3)
lines(c(25, 40, 100, 200, 300, 350, 450), em_preds1[1:7])
lines(c(25, 40, 100, 200, 300, 350, 450), em_preds1[8:14])
lines(c(25, 40, 100, 200, 300, 350, 450), em_preds1[1:7]-3*sqrt(em_vars1[1:7]), lty = 2)
lines(c(25, 40, 100, 200, 300, 350, 450), em_preds1[1:7]+3*sqrt(em_vars1[1:7]), lty = 2)
lines(c(25, 40, 100, 200, 300, 350, 450), em_preds1[8:14]-3*sqrt(em_vars1[8:14]), lty = 2)
lines(c(25, 40, 100, 200, 300, 350, 450), em_preds1[8:14]+3*sqrt(em_vars1[8:14]), lty = 2)
lines(c(25, 40, 100, 200, 300, 350, 450), em_preds2[1:7], col = 'red')
lines(c(25, 40, 100, 200, 300, 350, 450), em_preds2[8:14], col = 'red')
lines(c(25, 40, 100, 200, 300, 350, 450), em_preds2[1:7]-3*sqrt(em_vars1[1:7]), lty = 2, col = 'red')
lines(c(25, 40, 100, 200, 300, 350, 450), em_preds2[1:7]+3*sqrt(em_vars1[1:7]), lty = 2, col = 'red')
lines(c(25, 40, 100, 200, 300, 350, 450), em_preds2[8:14]-3*sqrt(em_vars1[8:14]), lty = 2, col = 'red')
lines(c(25, 40, 100, 200, 300, 350, 450), em_preds2[8:14]+3*sqrt(em_vars1[8:14]), lty = 2, col = 'red')
legend('topleft', legend = c('Recovered', "Infected"), lty = 1, col = c(3,4), inset = c(0.05, 0.05))
```

Last but not least: implausibility. nth_implausible deals with the bimodality automatically.
The way it works is to move output-wise: it calculates the implausibility for an output for
each mode, and takes the MINIMUM across the modes (since we only reject a point if it is deemed
                                                   implausible in both modes). We then take this collection of minimised implausibilities and
maximise/nth-maximise.

```{r, fig.height=7, fig.width=9}
emulator_plot(subset_emulators(bimodal_emulators, output_names[-c(7,14)]), plot_type = 'nimp', targets = targets, params = c('alpha', 'sigma'))
```

Let's generate some new points, see what we get.
Main thing to note is that, from the point of view of the user, nothing changes in terms of function use.

```{r}
new_points <- generate_new_runs(subset_emulators(bimodal_emulators, output_names[-c(7,14)]), 200, targets)
```
