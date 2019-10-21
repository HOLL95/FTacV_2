fig, ax=plt.subplots(harm_class.num_harmonics, 1)
line_elements=[]
for i in range(0, harm_class.num_harmonics):
    l1,=(ax[i].plot(time_results, abs(test_harmonics[i,:]), lw=2))
    (ax[i].plot(time_results, abs(data_harmonics[i,:])))
    line_elements.append(l1)


axcolor = 'lightgoldenrodyellow'
slider_ax=[]
slider_ax_element=[]
for i in range(0, len(ramp_fit.optim_list)):
    slider_ax.append(plt.axes([0.25, 0.0+(i*0.02), 0.65, 0.01], facecolor=axcolor))
    slider_ax_element.append(Slider(slider_ax[i], ramp_fit.optim_list[i], param_bounds[ramp_fit.optim_list[i]][0], param_bounds[ramp_fit.optim_list[i]][1],ramp_means[i]) )



def update(val):
    params=np.zeros(len(ramp_fit.optim_list))
    for i in range(0, len(ramp_fit.optim_list)):
        params[i]=slider_ax_element[i].val
    test=ramp_fit.test_vals(params, likelihood="timeseries", test=False)
    test_harmonics=harm_class.generate_harmonics(time_results, test)
    for i in range(0, harm_class.num_harmonics):
        line_elements[i].set_ydata(abs(test_harmonics[i,:]))
    fig.canvas.draw_idle()

for i in range(0, len(ramp_fit.optim_list)):
    slider_ax_element[i].on_changed(update)


resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')







plt.show()
ramp_fit.bounds_val=1e5
print(means[0:3])

plt.plot(test*-1)
plt.plot(current_results)
plt.show()
