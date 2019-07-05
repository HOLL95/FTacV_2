#plt.plot(time_results, current_results, alpha=0.5)
l, = plt.plot(time_results, test, lw=2)
plt.plot(time_results, current_results, alpha=0.5)


axcolor = 'lightgoldenrodyellow'
slider_ax=[]
slider_ax_element=[]
for i in range(0, len(ramp_fit.optim_list)):
    slider_ax.append(plt.axes([0.25, 0.0+(i*0.02), 0.65, 0.01], facecolor=axcolor))
    slider_ax_element.append(Slider(slider_ax[i], ramp_fit.optim_list[i], param_bounds[ramp_fit.optim_list[i]][0], param_bounds[ramp_fit.optim_list[i]][1], means[i]))



def update(val):
    params=np.zeros(len(ramp_fit.optim_list))
    for i in range(0, len(ramp_fit.optim_list)):
        params[i]=slider_ax_element[i].val
    cdl_idx=ramp_fit.optim_list.index("E_0")
    print params[cdl_idx]
    l.set_ydata(ramp_fit.simulate(params,frequencies, "yes", "timeseries", "no" ))
    fig.canvas.draw_idle()

for i in range(0, len(ramp_fit.optim_list)):
    slider_ax_element[i].on_changed(update)


resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')







plt.show()
ramp_fit.bounds_val=1e5
print means[0:3]

plt.plot(test*-1)
plt.plot(current_results)
plt.show()

plt.plot(current_results)
plt.plot(np.array(test)*-1)
plt.show()
