from pylab import *
import seaborn as sns
import pandas as pd

import ringdown as rd

line_color = sns.color_palette()[4]

# make plots fit the LaTex column size but rescale them for ease of display
scale_factor = 2

# Get columnsize from LaTeX using \showthe\columnwidth
fig_width_pt = scale_factor*246.0
# Convert pts to inches
inches_per_pt = 1.0/72.27               
# Golden ratio
fig_ratio = (np.sqrt(5)-1.0)/2.0
fig_width = fig_width_pt*inches_per_pt
fig_height =fig_width*fig_ratio

figsize_column = (fig_width, fig_height)
figsize_square = (fig_width, fig_width)

fig_width_page = scale_factor*inches_per_pt*508.87
figsize_page = (fig_width_page, fig_height)

rcParams = {'figure.figsize': figsize_column}

# LaTex text font sizse in points (rescaled as above)
fs = scale_factor*9
fs_label = 0.8*fs


def plot_pair_violin(dfs, keys=None, sr=2048, figsize=(10,3),
                    tlim=None):
    figure(figsize=figsize)

    if keys is None:
        keys = dfs.keys()
        
    try:
        sr[0]
    except TypeError:
        sr = [sr, sr]

    df = pd.DataFrame()
    for k, sr in zip(keys, sr):
        rdf = dfs[k][dfs[k]['srate']==sr].copy()
        if tlim is not None:
            rdf = rdf[(rdf['$t_0/M$'] >= tlim[0]) & (rdf['$t_0/M$'] <= tlim[1])]
        rdf['run'] = "{} ({:.0f})".format(k, sr)
        df = df.append(rdf)
        
    split = len(df['run'].unique()) > 1

    g = sns.violinplot(x='$t_0/M$', y="$A_1$", hue="run",
                       data=df, palette="Set2", split=split,
                       inner="quartile")
    axhline(0, ls='--', c='k', alpha=0.5)
    new_xlabels = ['{:.2f}'.format(float(i.get_text())) for i in g.get_xticklabels()]
    xticks(g.get_xticks(), new_xlabels);
    legend(loc='upper left');
    return g

def plot_sigmas(dfs, keys=None, sr=2048, figsize=(10,3), ax=None,
                m_ref=69, chi_ref=0.69, trend=False, **kws):
    tM = m_ref*rd.qnms.T_MSUN
    _, tau = rd.qnms.get_ftau(m_ref, chi_ref, 1)
    tauM = tau / tM
    
    lkws = dict(fmt='.', capsize=4, lw=2, capthick=2, alpha=0.7)
    lkws.update(kws)

    if keys is None:
        keys = dfs.keys()

    if ax is None:
        fig, ax = subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
    for i, key in enumerate(keys):
        df = dfs[key].get(sr, pd.DataFrame())
        lab = lkws.pop('label', key)
        a = lkws.get('alpha', 1)
        if not df.empty:
            yerr = (df['med'] - df['lo'], df['hi'] - df['med'])
            l = ax.errorbar(df.index, df['med'], yerr=yerr, label=lab, **lkws)
            c = l.get_children()[0].get_color()
            # plot expo trendline
            if trend:
                A0 = df['med'][min(abs(df.index[df.index>=0]))]
                t = df.index.values
                ax.plot(t, A0*exp(-t/tauM), c=c, alpha=a/2, ls='--')
                ax.plot(t[t>=0], A0*exp(-t[t>=0]/tauM), c=c, alpha=a)
    ax.axhline(0, ls='--', c='k')
    ax.legend();
    ax.set_ylabel('$A_1$')
    ax.set_xlabel('$t_0/t_M$')
    ax.set_title('analysis sampling rate {}'.format(sr));
    return fig

def plot_amps(fit=None, x=None, y=None, truth=None, d=1, g=None, levels=[0.9, 0.5, 0.1],
              points=True, truth_kws=None, xlim=(0, 8E-21), ylim=(0, 14E-21), **kws):
    with sns.plotting_context('paper', font_scale=1.5):
        g = g or sns.JointGrid(x=[], y=[], xlim=xlim, ylim=ylim)
        g.x = 2*fit.posterior.A[:,::d,0].values.flatten() if x is None else x
        g.y = 2*fit.posterior.A[:,::d,1].values.flatten() if y is None else y
        # set style
        lws = kws.get('lws', linspace(1, 2, len(levels)))
        lkws = kws.get('lkws', dict(lw=lws[-1], ls=kws.get('ls', '-')))
        l, = g.ax_joint.plot([], [], label=kws.get('label', None),
                             c=kws.pop('c', kws.pop('color', None)),
                             **lkws)
        c = l.get_color()
        # plot
        if points:
            g.plot_joint(scatter, color=c, alpha=0.03, marker='.')
        g.plot_joint(rd.kdeplot_2d_clevels, colors=[c,], cmap=None, levels=levels,
                    linewidths=lws, linestyles=kws.get('ls', '-'), **kws)
        calpha = matplotlib.colors.to_rgba(c, kws.get('alpha', None))
        g.plot_marginals(sns.kdeplot, c=calpha, **lkws)
        
        if truth:
            tkws = dict(c=c, ls='--')
            tkws.update(truth_kws or {})
            if isinstance(truth, dict):
                a0 = truth.get('A', truth['a'])[0]
                a1 = truth.get('A', truth['a'])[1]
            else:
                a0, a1 = truth[:2]
            g.ax_joint.axvline(a0, **tkws)
            g.ax_joint.axhline(a1, **tkws)
            g.ax_joint.plot(a0, a1, marker='+', ms=10, mew=1.5, **tkws)
        g.set_axis_labels(r'$A_0$', r'$A_1$');
    return g

def plot_mchi(fit=None, x=None, y=None, truth=None, d=1, g=None, 
              levels=[0.9, 0.5, 0.1], points=True, truth_kws=None,
              xlim=(50, 100), ylim=(0, 1), marginals=True, **kws):
    with sns.plotting_context('paper', font_scale=1.5):
        g = g or sns.JointGrid(x=[], y=[], xlim=xlim, ylim=ylim)
        g.x = fit.posterior.M[:,::d].values.flatten() if x is None else x
        g.y = fit.posterior.chi[:,::d].values.flatten() if y is None else y
        # set style
        lws = kws.get('lws', linspace(1, 2, len(levels)))
        lkws = kws.get('lkws', dict(lw=lws[-1], ls=kws.get('ls', '-')))
        l, = g.ax_joint.plot([], [], label=kws.get('label', None),
                             c=kws.pop('c', kws.pop('color', None)),
                             **lkws)
        c = l.get_color()
        # plot
        if points:
            g.plot_joint(scatter, color=c, alpha=0.03, marker='.')
        g.plot_joint(rd.kdeplot_2d_clevels, colors=[c,], cmap=None, levels=levels,
                    linewidths=lws, linestyles=kws.get('ls', '-'), **kws)
        calpha = matplotlib.colors.to_rgba(c, kws.get('alpha', None))
        if marginals:
            g.plot_marginals(sns.kdeplot, color=calpha, alpha=kws.get('alpha', None), **lkws)
        
        if truth:
            tkws = dict(c=c, ls='--')
            tkws.update(truth_kws or {})
            
            M = truth.get('M', truth['m'])
            chi = truth.get('chi')
            g.ax_joint.axvline(M, **tkws)
            g.ax_joint.axhline(chi, **tkws)
            g.ax_joint.plot(M, chi, marker='+', ms=10, mew=1.5, **tkws)
        g.set_axis_labels(r'$M/M_\odot$', r'$\chi$');
    return g

def plot_dfdtau(fit, truth=None, d=1, g=None, levels=[0.9, 0.5, 0.1], points=True,
                  truth_kws=None, xlim=(-0.5, 0.5), ylim=(-0.5, 0.5), npoints=1000, **kws):
    with sns.plotting_context('paper', font_scale=1.5):
        g = g or sns.JointGrid(x=[], y=[], xlim=xlim, ylim=ylim)
        n = fit.posterior.df.shape[0]*fit.posterior.df.shape[1]
        ixs = np.random.choice(n, min(n, npoints))
        g.x = fit.posterior.df[:,::d].values.flatten()[ixs]
        g.y = fit.posterior.dtau[:,::d].values.flatten()[ixs]
        # set style
        lws = kws.get('lws', linspace(1, 2, len(levels)))
        lkws = kws.get('lkws', dict(lw=lws[-1], ls=kws.get('ls', '-')))
        l, = g.ax_joint.plot([], [], label=kws.get('label', None),
                             c=kws.pop('c', kws.pop('color', None)),
                             **lkws)
        c = l.get_color()
        # plot
        if points:
            g.plot_joint(scatter, color=c, alpha=0.03, marker='.')
        g.plot_joint(rd.kdeplot_2d_clevels, colors=[c,], cmap=None, levels=levels,
                     linewidths=lws, linestyles=kws.get('ls', '-'), **kws)
        calpha = matplotlib.colors.to_rgba(c, kws.get('alpha', None))
        g.plot_marginals(sns.kdeplot, c=calpha, **lkws)
        
        if truth:
            tkws = dict(c=c, ls='--')
            tkws.update(truth_kws or {})
            
            g.ax_joint.axvline(truth.get('M', truth['m']), **tkws)
            g.ax_joint.axhline(truth.get('chi', truth['chi']), **tkws)
            g.ax_joint.plot(truth['M'], truth['chi'], marker='+', markersize=10,
                            markeredgewidth=1.5, **tkws)
        g.set_axis_labels(r'$\delta f_1$', r'$\delta \tau_1$');
    return g
