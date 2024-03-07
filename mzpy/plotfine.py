# 长远目标：常用生信图汇集于此，并统一样式

from itertools import chain
# from matplotlib import cm
from matplotlib import colormaps as cm
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
from plotnine import *
from sklearn.decomposition import PCA


class PlotFine():
    # 绘图的统一接口，只需使用该类即可
    def __init__(self, base_theme = None,
                fontsize = 12,
                figure_size = (5, 5),
                dpi = 120):
        if base_theme == None:
            base_theme = theme_matplotlib()
        self.theme = (base_theme + 
            theme(
                axis_title         = element_text(size=fontsize*1.3),
                axis_text          = element_text(size=fontsize),
                legend_title       = element_text(size=fontsize*1.3),
                legend_text        = element_text(size=fontsize),
                legend_background  = element_blank(),
                legend_position    = 'right',
                figure_size        = figure_size,
                dpi                = dpi))
        self.fontsize = fontsize
        self.figure_size = figure_size
        self.dpi = dpi
    
    def bubble(self, df:pd.DataFrame,
               x:str       ='impact',
               y:str       = '_pFDR_',
               fill:str    = 'category',
               size:str    = '_n_match_',
               n_top:int   = 20,
               palette:str = 'Set1',
               save_to:str = None):
        if df.shape[0] > n_top:
            df = df.head(n_top)
        n_duplicated = df.duplicated(subset=[x,y]).value_counts().get(True, 0)
        if n_duplicated > 0:
            print(f'There are {n_duplicated} data points overlapping!')

        plot = (ggplot(df, aes(x=x, y=y, size=size, fill=fill)) +
                geom_point(alpha=0.45) +
                scale_fill_brewer(type='qualitative', palette=palette) +
                self.theme
        )
        if save_to:
            plot.save(save_to, transparent=True)
        return plot
    
    def lollipop(self, df:pd.DataFrame, x:str, y:str, fill:str=None,
                    palette:str = 'Set1', save_to:str = None):
        df = df.sort_values(by = x, ascending = True)
        # 棒棒糖图的纵坐标必须转换为因子，否则绘图不排序
        df[y] = pd.Categorical(df[y],
                    categories = df[y].unique(),
                    ordered = True)
        if fill:
            p = ggplot(df, aes(x, y, fill=fill))
        else:
            p = ggplot(df, aes(x, y))
        plot = (p+
                    geom_segment(aes(x=0, xend=x, y=y, yend=y))+
                    geom_point(shape='o', size=3, color='black')+
                    scale_fill_brewer(type='qualitative', palette=palette)+
                    self.theme)
        if save_to is not None:
            p.save(save_to, transparent=True)
        return plot

    def pca(self, data, groups:list = None, labels:list = None,
            palette='Set1',
            save_to:str = None):
        pca = PCA(n_components=2).fit(data)        
        df = pd.DataFrame(pca.transform(data), columns=['PC1', 'PC2'])
        df['group'] = pd.Categorical(groups)
        plot = (ggplot(df, aes('PC1', 'PC2', fill='group'))+
                geom_point(alpha = 0.6, size = 3, shape = 'o', stroke = 0)+            
                stat_ellipse(geom="polygon", level=0.95, alpha=0.2)+
                labs(x = "PC1: %.1f %%"%(100*pca.explained_variance_ratio_[0]),
                    y = "PC2: %.1f %%"%(100*pca.explained_variance_ratio_[1]))+
                scale_fill_brewer(type='qualitative', palette=palette)+
                self.theme
        )
        if labels is not None:
            df['label'] = labels
            plot = plot + geom_text(label=df.label, nudge_x=0.1, nudge_y=0.1,
                              size = self.fontsize*0.6)
        if save_to:
            plot.save(save_to, transparent=True)
        return plot

    def volcano(self,df, x, y, fill, xcut = 1, ycut = 2,
                title = '',
                palette='Set1', save_to=None):
        x_limit = max(abs(min(df[x])), abs(max(df[x]))) 
        plot = (ggplot(df, aes(x, y, fill=fill)) +
                geom_point (alpha=0.5, size=2, shape='o', stroke=0) +
                geom_vline(xintercept =  xcut, linetype = "dashed", color = "gray", alpha = 0.7) +
                geom_vline(xintercept = -xcut, linetype = "dashed", color = "gray", alpha = 0.7) +            
                geom_hline(yintercept =  ycut, linetype = "dashed", color = "gray", alpha = 0.7) +
                xlim(-x_limit, x_limit) +
                ggtitle(title) +
                scale_fill_brewer(type='qualitative', palette=palette)+
                self.theme
        )
        if save_to:
            plot.save(save_to, transparent=True)
        return plot

    def venn(self, data, palette='Set1', alpha=0.65, save_to:str = None):
        plot = Venn(data,
                    palette = palette,
                    fontsize = self.fontsize,
                    alpha=alpha,
                    save_to = save_to)
        return plot



class Venn(): # 类名首字母是大写
    ''' plot venn
        modified on pyvenn: https://github.com/tctianchi/pyvenn
        kepping venn2, venn3 and venn4
    '''    
    def __init__(self, data, fill=['number'],
                 palette='Set1', fontsize=14, alpha=0.65,
                 save_to=None):
        # finished venn plot just here
        plt.clf() # 清除当前图形
        self.fig = plt.figure(0, figsize=(9,7), dpi=96)
        self.ax  = self.fig.add_subplot(111, aspect='equal')
        self.ax.set_axis_off()
        self.ax.set_ylim(bottom=0.0, top=1.0)
        self.ax.set_xlim(left=0.0, right=1.0)

        self.palette = palette
        self.fontsize = fontsize
        self.alpha = alpha
        n = len(data)
        values = list(data.values())
        keys   = list(data.keys())
        labels = self.divide(values, fill=fill)
        if n in (2, 3, 4):
            func = eval('self._venn%d'%n)
        else:
            raise ValueError('The data length must be in 2, 3 or 4')
        func(labels, keys) ### 统一接口
        if save_to:
            # plt.savefig(save_to, transparent=True)
            self.fig.savefig(save_to, transparent=True)

    def divide(self, data, fill=["number"]):
        """get a dict of labels for groups in data
        @type data: list[Iterable]
        @rtype: dict[str, str]
        input
            data: data to get label for
            fill: ["number"|"logic"|"percent"]
        return
            labels: a dict of labels for different sets
        example:
        In [12]: ([range(10), range(5,15), range(3,8)], fill=["number"])
        Out[12]:
        {'001': '0', '010': '5', '011': '0', '100': '3', '101': '2',
        '110': '2',  '111': '3'}
        """
        N = len(data)
        sets_data = [set(data[i]) for i in range(N)] # sets for separate groups
        s_all = set(chain(*data))                    # union of all sets
        set_collections = {}
        for n in range(1, 2**N):
            key = bin(n).split('0b')[-1].zfill(N)
            value = s_all
            sets_for_intersection = [sets_data[i] for i in range(N) if key[i] == '1']
            sets_for_difference = [sets_data[i] for i in range(N) if key[i] == '0']
            for s in sets_for_intersection:
                value = value & s
            for s in sets_for_difference:
                value = value - s
            set_collections[key] = value
        labels = {k: "" for k in set_collections}
        if "logic" in fill:
            for k in set_collections:
                labels[k] = k + ": "
        if "number" in fill:
            for k in set_collections:
                labels[k] += str(len(set_collections[k]))
        if "percent" in fill:
            data_size = len(s_all)
            for k in set_collections:
                labels[k] += "(%.1f%%)" % (100.0 * len(set_collections[k]) / data_size)
        return labels

    def _venn2(self, labels, names=['A', 'B']):
        """
        plots a 2-set Venn diagram

        @type labels: dict[str, str]
        @type names: list[str]
        @rtype: (Figure, AxesSubplot)

        input
            labels: a label dict where keys are identified via binary codes ('01', '10', '11'),
                    hence a valid set could look like: {'01': 'text 1', '10': 'text 2', '11': 'text 3'}.
                    unmentioned codes are considered as ''.
            names:    group names
            more:     colors, figsize, dpi, fontsize

        return
            pyplot Figure and AxesSubplot object
        """
        self.ax.set_ylim(bottom=0.0, top=0.7)
        # body
        xy = [((0.375, 0.3)), (0.625, 0.3)]
        for i in range(2):
            self.ax.add_patch(patches.Ellipse(xy     = xy[i],
                                        width  = 0.5,
                                        height = 0.5,
                                        alpha  = self.alpha,
                                        color  = cm[self.palette](i)))
        text = [(0.74, 0.30, labels.get('01', '')),
                (0.26, 0.30, labels.get('10', '')), 
                (0.50, 0.30, labels.get('11', '')),
                (0.20, 0.56, names[0]),
                (0.80, 0.56, names[1])]
        for it in text:
            self.ax.text(it[0], it[1],it[2], fontsize=self.fontsize, ha='center', va='center')

    def _venn3(self, labels, names=['A', 'B', 'C']):
        # body
        xy = [(0.333, 0.633), (0.666, 0.633), (0.500, 0.310)]
        for i in range(3):
            self.ax.add_patch(patches.Ellipse(xy = xy[i],
                                        width  = 0.55,
                                        height = 0.55,
                                        alpha  = self.alpha,
                                        color  = cm[self.palette](i)))
        text = [(0.50, 0.27, labels.get('001', '')),
                (0.73, 0.65, labels.get('010', '')),
                (0.61, 0.46, labels.get('011', '')),
                (0.27, 0.65, labels.get('100', '')),
                (0.39, 0.46, labels.get('101', '')),
                (0.50, 0.65, labels.get('110', '')),
                (0.50, 0.51, labels.get('111', '')),
                (0.15, 0.87, names[0]),
                (0.85, 0.87, names[1]),
                (0.50, 0.02, names[2])]
        for it in text:
            self.ax.text(it[0], it[1],it[2], fontsize=self.fontsize, ha='center', va='center')

    def _venn4(self, labels, names=['A', 'B', 'C', 'D']):
        o  = 0.500 # 图形中心点坐标：x=y=0.5
        dx = 0.18
        dy = 0.08
        xy = [(o-dx, o-dy), (o, o), (o, o), (o+dx, o-dy)]
        for i in range(4):
            angle = 135 if i < 2 else 45
            self.ax.add_patch(patches.Ellipse(xy     = xy[i],
                                        width  = 4*dx,
                                        height = 2*dx,
                                        angle  = angle,
                                        alpha  = self.alpha,
                                        color  = cm[self.palette](i)))
        text = [(o+dx*2.00, o+dy*0.50, labels.get('0001', '')),
                (o+dx*0.75, o+dy*2.50, labels.get('0010', '')),
                (o+dx*1.25, o+dy*1.25, labels.get('0011', '')),
                (o-dx*0.75, o+dy*2.50, labels.get('0100', '')),
                (o+dx     , o-dy*2.00, labels.get('0101', '')),
                (o        , o+dy*1.25, labels.get('0110', '')),
                (o+dx*0.75, o-dy*0.25, labels.get('0111', '')),
                (o-dx*2.00, o+dy*0.50, labels.get('1000', '')),
                (o        , o-dy*3.75, labels.get('1001', '')),
                (o-dx     , o-dy*2.00, labels.get('1010', '')),
                (o-dx*0.25, o-dy*2.75, labels.get('1011', '')),
                (o-dx*1.25, o+dy*1.25, labels.get('1100', '')),
                (o+dx*0.25, o-dy*2.75, labels.get('1101', '')),
                (o-dx*0.75, o-dy*0.25, labels.get('1110', '')),
                (o        , o-dy*1.75, labels.get('1111', '')),
                (o-dx*2.25, o+dy*2.75, names[0]),
                (o-dx*1.00, o+dy*3.75, names[1]),
                (o+dx*1.00, o+dy*3.75, names[2]),
                (o+dx*2.25, o+dy*2.75, names[3])]
        for it in text:
            self.ax.text(it[0], it[1], it[2], fontsize=self.fontsize, ha='center', va='center')

