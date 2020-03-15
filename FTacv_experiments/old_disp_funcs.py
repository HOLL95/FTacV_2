else:
            if ("E0_std" in self.optim_list) and ("k0_shape" in self.optim_list):
                e0_vals, e0_disp=self.therm_dispersion()
                k0_vals, k0_disp=self.kinetic_dispersion()
                values=list(itertools.product(e0_vals, k0_vals))
                flags=list(zip(["E_0"]*len(values), ["k_0"]*len(values)))
                weights=list(itertools.product(e0_disp, k0_disp))
                weights=[weights[i][0]*weights[i][1] for i in range(len(weights))]
                weight_val_tuple=list(zip(flags, values, weights))
                paralell=paralell_class(self.nd_param_dict, self.time_vec, "sinusoidal", self.bounds_val, isolver_martin_brent.brent_current_solver)
                time_series=paralell.paralell_dispersion(weight_val_tuple)

            elif "E0_std" in self.optim_list:
                start1=time.time()
                e0_vals, e0_disp=self.therm_dispersion()

                counter=0
                if "alpha_dispersion" in self.simulation_options:
                    if self.simulation_options["alpha_dispersion"]=="uniform":
                        alpha_bins=5
                        alpha_vals=np.linspace(self.param_bounds["alpha"][0],self.param_bounds["alpha"][1], alpha_bins)
                        alpha_weights=[1/alpha_bins]*alpha_bins
                    elif self.simulation_options["alpha_dispersion"]=="normal":
                        if "GH_quadrature" in self.simulation_options:
                            if self.simulation_options["GH_quadrature"]==True:
                                alpha_vals=self.normal_gh_transform(location=self.nd_param.alpha_mean, scale=self.nd_param.alpha_std)
                                alpha_weights=self.normal_GH_weights
                                alpha_bins=len(alpha_vals)
                        else:
                            alpha_bins=16
                            alpha_min=max(0, norm.ppf(1e-4, loc=self.nd_param.alpha_mean, scale=self.nd_param.alpha_std))
                            alpha_max=min(1,norm.ppf(1-(1e-4), loc=self.nd_param.alpha_mean, scale=self.nd_param.alpha_std))
                            alpha_weights=np.zeros(alpha_bins)
                            alpha_vals=np.linspace(alpha_min, alpha_max, alpha_bins)
                            alpha_weights[0]=norm.cdf(alpha_vals[0], loc=self.nd_param.alpha_mean, scale=self.nd_param.alpha_std)
                            alpha_mids=np.zeros(alpha_bins)
                            alpha_mids[0]=alpha_vals[0]
                            for i in range(1, len(alpha_weights)):
                                alpha_weights[i]=norm.cdf(alpha_vals[i],loc=self.nd_param.alpha_mean, scale=self.nd_param.alpha_std)-norm.cdf(alpha_vals[i-1],loc=self.nd_param.alpha_mean, scale=self.nd_param.alpha_std)
                                alpha_mids[i]=(alpha_vals[i]+alpha_vals[i-1])/2
                            alpha_vals=alpha_mids
                    for i in range(0, alpha_bins):
                        self.nd_param_dict["alpha"]=float(alpha_vals[i])
                        for j in range(0,self.simulation_options["dispersion_bins"]):
                            self.nd_param_dict["E_0"]=float(e0_vals[j])
                            time_series_current=solver(self.nd_param_dict, self.time_vec,self.simulation_options["method"], -1, self.bounds_val)
                            time_series=np.add(time_series, np.multiply(time_series_current, e0_disp[j]*alpha_weights[i]))
                else:
                    for i in range(0, self.simulation_options["dispersion_bins"]):
                        self.nd_param_dict["E_0"]=float(e0_vals[i])
                        time_series_current=solver(self.nd_param_dict, self.time_vec,self.simulation_options["method"], -1, self.bounds_val)
                        time_series=np.add(time_series, np.multiply(time_series_current, e0_disp[i]))
                        #plt.plot(voltages, np.multiply(time_series_current, e0_disp[i]))
                        #print(e0_disp[i])
            elif ("k0_shape" in self.optim_list):
                k0_vals, k0_disp=self.kinetic_dispersion()
                for i in range(0, self.simulation_options["dispersion_bins"]):
                    self.nd_param_dict["k_0"]=k0_vals[i]
                    time_series_current=solver(self.nd_param_dict, self.time_vec,self.simulation_options["method"], -1, self.bounds_val)
                    time_series=np.add(time_series, np.multiply(time_series_current, k0_disp[i]))
            #plt.plot(k0_vals, k0_disp)
            #plt.show()

            def kinetic_dispersion(self):
        #print self.nd_param.k0_shape, self.nd_param.k0_loc, self.nd_param.k0_scale
        k0_weights=np.zeros(self.simulation_options["dispersion_bins"])
        k_start=lognorm.ppf(1e-9, self.nd_param.k0_shape, 0, scale=self.nd_param.k0_scale)
        k_end=lognorm.ppf(1-1e-5, self.nd_param.k0_shape, 0, scale=self.nd_param.k0_scale)
        k0_vals=np.linspace(k_start,k_end, self.simulation_options["dispersion_bins"])
        k0_weights[0]=lognorm.cdf(k0_vals[0], self.nd_param.k0_shape, 0, scale=self.nd_param.k0_scale)
        for k in range(1, len(k0_weights)):
            k0_weights[k]=lognorm.cdf(k0_vals[k], self.nd_param.k0_shape, 0, scale=self.nd_param.k0_scale)-lognorm.cdf(k0_vals[k-1], self.nd_param.k0_shape, 0, scale=self.nd_param.k0_scale)
        #plt.title("k0")
        return k0_vals, k0_weights

    def therm_dispersion(self):
        if "GH_quadrature" in self.simulation_options:
            if self.simulation_options["GH_quadrature"]==True:
                e0_vals=self.normal_gh_transform(location=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
                e0_weights=self.normal_GH_weights

        else:
            self.e0_min=norm.ppf(1e-4, loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
            self.e0_max=norm.ppf(1-(1e-4), loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
            e0_weights=np.zeros(self.simulation_options["dispersion_bins"])
            e0_vals=np.linspace(self.e0_min, self.e0_max, self.simulation_options["dispersion_bins"])
            e0_midpoints=np.zeros(self.simulation_options["dispersion_bins"])
            e0_weights[0]=norm.cdf(e0_vals[0], loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
            e0_midpoints[0]=(e0_vals[0]+self.e0_min)/2
            for i in range(1, len(e0_weights)):
                e0_weights[i]=norm.cdf(e0_vals[i],loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)-norm.cdf(e0_vals[i-1],loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
                e0_midpoints[i]=(e0_vals[i]+e0_vals[i-1])/2
            e0_vals=e0_midpoints
            #midpoint=len(e0_weights)//2
            #e0_weights[midpoint:]=e0_weights[:midpoint]
            #range1=np.arange(self.simulation_options["dispersion_bins"], 0, -1)
            #e0_weights=np.divide(range1, sum(range1))
            #print(e0_weights)
            #plt.plot(self.e_nondim(e0_vals), e0_weights)
            #plt.show()
        #for i in range(0, len(e0_midpoints)):
        #    plt.axvline(self.e_nondim(e0_midpoints[i]), color="black", linestyle="--")
        #plt.show()

        #print(list(e0_weights))
        #print(list(self.e_nondim(e0_vals)))
        #plt.axvline(self.e_nondim(self.nd_param.E0_mean))
        #plt.plot(self.e_nondim(e0_vals), e0_weights)
        #print self.nd_param.E0_mean,self.nd_param.E0_std
        #plt.title("e0")
        #plt.show()
        return e0_vals, e0_weights
    def alpha_dispersion(self):
        self.e0_min=norm.ppf(1e-11, loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
        self.e0_max=norm.ppf(1-(1e-11), loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
        e0_weights=np.zeros(self.simulation_options["dispersion_bins"])
        e0_vals=np.linspace(self.e0_min, self.e0_max, self.simulation_options["dispersion_bins"])
        e0_weights[0]=norm.cdf(e0_vals[0], loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
        for i in range(1, len(e0_weights)):
            e0_weights[i]=norm.cdf(e0_vals[i],loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)-norm.cdf(e0_vals[i-1],loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
        #plt.plot(e0_vals, e0_weights)
        #print self.nd_param.E0_mean,self.nd_param.E0_std
        #plt.title("e0")
        #plt.show()
        return e0_vals, e0_weights
    def weight_matrix(self,e0_disp, k0_disp):
        e0_mat, k0_mat=np.meshgrid(e0_disp, k0_disp)
        weights=np.multiply(e0_mat, k0_mat)
        return weights
