# Saving Results
GridSearchSys objects can be saved to file with the Serialization library and a bit of overhead work, encapsulated in `save_serde_data` and `load_serde_data`.
```@docs
save_serde_data
load_serde_data
```
They both save to or load from a *directory*, since they work with multiple files.

Saving to file lets you do neat things like

```julia
gss = load_serde_data("path/to/gss/folder")
add_result!(gss, "eigenvalues", eu.get_eigs)
save_serde_data("path/to/gss/folder")
```

In addition, the `set_chunksize!` method controls chunking.

```@docs
set_chunksize!
```

The default chunk size is `Inf`. If this is the case, then when `save_serde_data` is called, it will just save the entire dataframe of results to a single file. Otherwise, it will chunk it into files with `gss.chunksize` rows.


When `execute_sims!` is called, it writes the results to `gss.df` one row at a time. When it reaches `gss.chunksize` total rows, it saves the entire dataframe to file, then deletes it from RAM. This allows very very large sweeps to be run in finite memory.

When it finishes all of the simulations, it saves whatever is left to file as well.

To prevent it from saving anything, you can run

```julia
set_chunksize(gss, Inf)
```

## The header file and associated shenanigans
One issue with serialization is that it doesn't work with user-defined functions. For example, if we tried to save our `gss` from before, it would save some reference to the `set_power_setpt!` function, but wouldn't save it. Therefore, we need to redefine it wherever we want to load back our results. 

But unless you're trying to run the simulations again, you don't *really* need the function definition. You just need something for it to refer to. As a result, my quick hack is that every time we reference an external function, we add `"function {name} end; "` to the `gss.hfile` variable. Then we eval this string when we load the `gss` object back, and it all works!

...In theory, at least. The best solution is always to simply import or copy the definitions into wherever you need them. If you are having issues, you can typically just open the file (`.hfile` in the save directory) and delete whichever functions are screaming at you.

## If all else fails
Obviously, the saved data is of great importance, since it often takes forever to produce. Therefore, I wanted to make sure that even if my funny additions on top of `Serialization` didn't work, you wouldn't lose your data.

Accordingly, the `load_serde_data` function will also accept a path to one of the `.jls` files saved by `save_serde_data`. These are just plain serialized `DataFrame`s, and they load back quite reliably. 

If you don't *have* to do this, avoid it, but it's there just in case.