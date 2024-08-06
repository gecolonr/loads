# Plotting

This one's complicated. Yikes.

Essentially, build your plots with `makeplots` and save them to interactive html with `savehtmlplot`.

For `makeplots`, you pass the dataframe (typically `gss.df` or similar), a bunch of Strings or Symbols specifying columns to map to various variables, and whatever other random settings you need.

The docstring is horrendously complicated and so is the signature, but remember: if you don't need the functionality, just don't include it, and it'll (maybe (probably (surely))) be fine.

```@docs
savehtmlplot
makeplots
```

## Things to be aware of
Also known as Questionable Things I'm Not Willing To Fix Right Now

 - if `colorbar==false`, `color` is assumed to be categorical and each unique color is mapped to an element of `colorlist` after being sorted by `color_sort_func`
 - colorbar doesn't really work with lines because you can't have lines of variable color. To get it to work, you need to pass `colorlist` and ensure that it matches the right values on the colorscale you give in `colorbar_args`
 - `colorbar_args` are kinda complicated. Here's an example:
```julia
Dict(attr(
    autocolorscale=false, 
    cmax= 1.0,
    # cmid=0.5,
    cmin= -1.2,
    title="",
    colorbar=attr(
        outlinecolor=colorant"black",
        outlinewidth=1,
        tickmode="array",
        ticktext=(["Z=1.0", "P=1.0", "P=0.8", "P=0.6", "P=0.4", "P=0.2", "P=E=0", "E=0.2", "E=0.4", "E=0.6", "E=0.8", "E=1.0"]),
        tickvals=(collect(range(-1.2, 1, 13)).+(2.2/(2*12)))[1:end - 1],
        xref="paper",
        yref="paper",
        x=1.02,
        thickness=40
    ),
    colorscale=[[0.0, RGB(colorant"yellow")], 
                [1.0, RGB(colorant"purple")]]
))
```
Note that `colorscale`'s numerical scale is in $$[0, 1]$$, which is different from `cmax` and `cmin`.
 - You can have a legend and a colorbar at the same time, but it is recommended to just use `hide_legend_duplicates=true` (which is default `true`) and use the legend to represent a *different* variable. Since color cannot be used, you'll have to pair it with another method of conveying information, like opacity or marker shape.
 - `symbollist` and `colorlist` provide some initial values, but if you have more unique values represented by symbol or color respectively, you'll need to provide new lists with more new values.
 - `scattertext` and `hovertext` can make plot files MASSIVE. Take care to include as little information as possible in them. Round numbers to the appropriate length.
 - `slider_trace_id` is funny. You don't actually have to make some fancy id - often using the default, the order of the rows in the dataframe, is best. Just make sure to sort the dataframe beforehand by as many things as possible so that you have a very predictable and consistent order of traces across slider values.