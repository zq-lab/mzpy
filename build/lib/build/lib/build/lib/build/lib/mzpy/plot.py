from itertools import chain

from matplotlib import colormaps as cm
from matplotlib.colors import ListedColormap, to_hex
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

import pandas as pd
from plotnine import *
import seaborn as sns
from sklearn.decomposition import PCA

class Plot():
    def __init__(self, base_theme = None,
                fontsize = 20,
                figure_size = (5, 5),
                dpi = 96):
        if base_theme == None:
            base_theme = theme_matplotlib()

        self.theme = (base_theme + 
            theme(
                axis_title         = element_text(size=fontsize*1.3),
                axis_text          = element_text(size=fontsize),
                legend_title       = element_text(size=fontsize),
                legend_text        = element_text(size=fontsize*0.85),
                legend_background  = element_rect(fill=None, color=None),
                legend_position    = (0.95, 0.95),
                figure_size        = figure_size,
                dpi                = dpi))
        
        self.fontsize = fontsize
        self.figure_size = figure_size
        self.dpi = dpi

        # Nature 系列参考 https://zhuanlan.zhihu.com/p/670396774
        self.colors = {
            'Nature_1': ['#217185', '#D95319', '#FED976', '#77AC30'],
            'Nature_2': ['#A5AEB7', '#925EB0', '#37E99F4', '#CC7C71', '#7AB656']
        }
        
    
    def bubble(self, df,
               x:str       ='impact',
               y:str       = '_pFDR_',
               fill:str    = 'category',
               size:str    = '_n_match_',
               n_top:int   = None,
               palette:str = 'Dark2',
               save_to:str = None):
        # 需要注意的是，像 'tab10' 这种来自 Matplotlib 的命名通常不在原生的 ColorBrewer 范围里
        # 可能会出现找不到该调色板或默认回退为其他 palette 的情况。
        if n_top:
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
    
    
    def decision_regions(self, X, y, labels, classifier, resolution=0.02, cmap='tab10', save_to=None):
        '''
        绘制决策区域
        X, feature matrix
        y, lable vector,correspondings to X
        classfier, for predict decisiong regions,in most cased, it would be an object of LogisticRegression
        labels, 
        resolution, meshgrid resolution
        cmap, cmap for dot of samples.
        '''
        markers = ('s', 'x', 'o', '^', 'v')
        colors = cm[cmap].colors

        # plot the decision surface
        x1_min, x1_max = X[:, 0].min() - 1, X[:, 0].max() + 1
        x2_min, x2_max = X[:, 1].min() - 1, X[:, 1].max() + 1
        xx1, xx2 = np.meshgrid(np.arange(x1_min, x1_max, resolution),
                            np.arange(x2_min, x2_max, resolution))
        Z = classifier.predict(np.array([xx1.ravel(), xx2.ravel()]).T)
        Z = Z.reshape(xx1.shape)
        plt.figure(figsize=(5, 5))
        plt.contourf(xx1, xx2, Z, alpha=0.4, cmap='Pastel2', antialiased=True)
        plt.xlim(xx1.min(), xx1.max())
        plt.ylim(xx2.min(), xx2.max())

        # plot samples by class
        for idx, cl in enumerate(np.unique(y)):
            plt.scatter(x=X[y == cl, 0], 
                        y=X[y == cl, 1],
                        alpha=0.8, 
                        color=colors[idx],
                        marker=markers[idx], 
                        label=labels[idx])
        plt.legend(loc='lower left')
        plt.tight_layout()
        if save_to:
            plt.savefig(save_to, transparent=True)
        plt.show()


    def heatmap(self, df,
                data_transfer = 'log10',
                palette='seismic', 
                title = r"Heatmap of Log$_{10}$ (Peak Area)",
                xlab = "Sample",
                ylab = "Metabolite",
                save_to = None):
        #  将数据转换为对数（自然对数），使用 apply 方法
        if data_transfer == 'log10':
            data = df.apply(lambda x: np.log(x + 1))  # 加1避免对0取对数
        elif data_transfer == 'relative':
            base_pk = df.max().max()  # 获取整个 DataFrame 的最大值  
            data = (df / base_pk) * 100  
        elif (data_transfer is None) or (data_transfer == ''):
            data = df
        else:
            raise ValueError(f'unknown data_transfer was set: {data_transfer}. (log10, relative or None are acceptable.)')

        # 设置绘图风格（移除网格线）
        sns.set_style("white")

        # 动态计算图像大小  
        num_rows, num_cols = data.shape  
        # 设置图像大小，确保小方块在不同数据量情况下保持近似大小  
        plt.figure(figsize=(max(5, num_cols * 0.5), max(5, num_rows * 0.5)))  # 根据列数和行数动态调整图像大小
        plot = sns.heatmap(
            data,
            cmap=palette,  # 使用蓝-白-红渐变色彩方案
            center=np.median(data.values.flatten()),  # 将中间值设置为中位数
            linewidths=0.05,
            linecolor='white',
            square=True,  # 单元格强制为正方形
            cbar_kws={"shrink": 0.8, "label": r'log$_{10}$ (peak area)'}  # 调整颜色条大小
        )

        # 设置标题和标签
        plt.title(title, fontsize=16)
        plt.xlabel(xlab, fontsize=12)
        plt.ylabel(ylab, fontsize=12)

        # # 调整 x 轴标签的旋转角度
        plot.set_xticklabels(plot.get_xticklabels(), rotation=90, ha="center")
        plot.set_yticklabels(plot.get_yticklabels(), rotation=0, ha="right")  # 设置为0度，右对齐

        if save_to:
            plt.savefig(save_to, format=save_to.split('.')[-1], bbox_inches='tight')    

        plt.show() 

    def line_with_error_band(self,
                             df_long,
                             id_on = 'id',
                             group_on = 'group',
                             value_on='peak area',
                             palette='tab10',
                             save_to=None):
        '''
        绘制带误差线的折线图
        param:
            df_long, long table containing id, group and values (value_name)
            by, column names list (or name) which to calculate mean and std
            value_on, column name of values (maybe peak area or height, or normalized values) 
        '''

        # 计算每个 group 和 id 的均值和标准差  
        summary = df_long.groupby([id_on, group_on]).agg(  
            Mean=(value_on, 'mean'),  
            StdDev=(value_on, 'std')  
        ).reset_index()  

        # 计算误差带  
        summary['Error'] = summary['StdDev']  

        cmap = plt.colormaps[palette]  # 获取调色板的前两个颜色  
        line_color = to_hex(cmap(0))  # 第一个颜色  
        ribbon_color = to_hex(cmap(1))  # 第二个颜色  

        # 绘制折线图  
        plot = (  
            ggplot(summary, aes(x=group_on, y='Mean', group=id_on)) +  # 去掉颜色映射  
            geom_line(color=line_color) +  # 设置线条颜色为调色板的第一个颜色  
            geom_ribbon(aes(ymin='Mean - Error', ymax='Mean + Error'), alpha=0.1, fill=ribbon_color) +  # 设置误差带颜色  
            facet_wrap(f'~ {id_on}', scales='free_y') +  
            labs(title=f'{value_on} by {id_on} with Error Bands', x=group_on, y=value_on) +  
            theme_bw() +  
            theme(legend_position='none',
                  aspect_ratio=1,
                  panel_grid_major=element_blank(),  # 去掉主要网格线  
                  panel_grid_minor=element_blank()   # 去掉次要网格线
                ) +  
            coord_fixed()  # 确保每个子图保持正方形  
        )  

        if save_to:
            plot.save(save_to, transparent=True)

        return plot 

    
    def lloly(self, df:pd.DataFrame, x:str, y:str, fill:str=None,
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
            plot.save(save_to, transparent=True)
        return plot
    

    # def pca(self, data, groups:list = None, labels:list = None,
    #         palette='Set1',
    #         save_to:str = None):
    #     pca = PCA(n_components=2).fit(data)        
    #     df = pd.DataFrame(pca.transform(data), columns=['PC1', 'PC2'])
    #     df['group'] = pd.Categorical(groups)
    #     plot = (ggplot(df, aes('PC1', 'PC2', fill='group'))+
    #             geom_point(alpha = 0.6, size = 3, shape = 'o', stroke = 0)+            
    #             stat_ellipse(geom="polygon", level=0.95, alpha=0.2)+
    #             labs(x = "PC1: %.1f %%"%(100*pca.explained_variance_ratio_[0]),
    #                 y = "PC2: %.1f %%"%(100*pca.explained_variance_ratio_[1]))+
    #             scale_fill_brewer(type='qualitative', palette=palette)+
    #             self.theme
    #     )
    #     if labels is not None:
    #         df['label'] = labels
    #         plot = plot + geom_text(label=df.label, nudge_x=0.1, nudge_y=0.1,
    #                           size = self.fontsize*0.6)
    #     if save_to:
    #         plot.save(save_to, transparent=True)
    #     return plot

    def pca(self, data, groups: list = None, labels: list = None,
            palette='Set2',
            draw_ellipse: bool = True,
            min_samples_for_ellipse: int = 4, # 如果样本数太少，比如3个样本，plotnine不会绘制椭圆
            save_to: str = None):
        '''
        data, a matrix, rows are samples, columns are features (genes, proteins, or metabolites)        
        '''

        pca = PCA(n_components=2).fit(data)        
        df = pd.DataFrame(pca.transform(data), columns=['PC1', 'PC2'])
        df['group'] = pd.Categorical(groups)
        
        plot = (ggplot(df, aes('PC1', 'PC2', fill='group'))+
                geom_point(alpha = 0.6, size = 3, shape = 'o', stroke = 0)+            
                labs(x = "PC1: %.1f %%"%(100*pca.explained_variance_ratio_[0]),
                    y = "PC2: %.1f %%"%(100*pca.explained_variance_ratio_[1]))+
                scale_fill_brewer(type='qualitative', palette=palette)+
                self.theme
        )
        
        # 检查是否绘制椭圆
        if draw_ellipse and groups is not None:
            # 检查每个分组的样本数量
            group_counts = pd.Series(groups).value_counts()
            min_count = group_counts.min()
            
            if min_count >= min_samples_for_ellipse:
                plot = plot + stat_ellipse(geom="polygon", level=0.95, alpha=0.2)
            else:
                info = f'Warning: Minimum number of samples per group is {min_count}, '
                info = info + f'less than {min_samples_for_ellipse}, skip drawing the ellipse'
                print(info)
        
        if labels is not None:
            df['label'] = labels
            plot = plot + geom_text(label=df.label, nudge_x=0.1, nudge_y=0.1,
                            size = self.fontsize*0.6)
        if save_to:
            plot.save(save_to, transparent=True)
        return plot
    
    def pca_sns(self, data, groups: list = None, labels: list = None,
                palette='Set2',
                draw_ellipse: bool = True,
                min_samples_for_ellipse: int = 3,
                save_to: str = None):
        '''
        data, a matrix, rows are samples, columns are features (genes, proteins, or metabolites)
        pca函数使用plotnine绘图，当样本数低于4时，无法绘制置信椭圆
        本函数采用sns绘图，手动添加上置信椭圆图，即使样本数低于4，也可绘制椭圆图
        
        '''
        
        from matplotlib.patches import Ellipse
        from scipy.stats import chi2

        # 1) PCA分析
        pca = PCA(n_components=2, random_state=42)
        scores = pca.fit_transform(data)  # 假设 data 是特征矩阵

        # 2) 组装 DataFrame 以便按组绘图
        df = pd.DataFrame(scores, columns=["PC1", "PC2"])
        if groups is not None:
            df["condition"] = groups
        else:
            df["condition"] = labels

        expl = pca.explained_variance_ratio_  # PC1和PC2的解释度
        print(f'Interpretability: {expl}')

        # 3) 置信椭圆函数
        def plot_confidence_ellipse(data2d, ax, edgecolor="black", facecolor=None,
                                    alpha=0.15, linewidth=1.0, z=0.95):
            if data2d.shape[0] < 2:
                return
            cov = np.cov(data2d, rowvar=False)
            mean = data2d.mean(axis=0)

            # 特征分解获取主轴与方差
            vals, vecs = np.linalg.eigh(cov)
            order = vals.argsort()[::-1]
            vals, vecs = vals[order], vecs[:, order]

            # 椭圆旋转角度（度）
            theta = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))

            # 2维卡方阈值：sqrt(chi2.ppf(z, df=2))
            chi2_val = chi2.ppf(z, df=2)
            width, height = 2 * np.sqrt(vals * chi2_val)  # 直径 = 2×半轴

            ell = Ellipse(
                xy=mean, width=width, height=height, angle=theta,
                edgecolor=edgecolor, facecolor=facecolor if facecolor is not None else edgecolor,
                alpha=alpha, linewidth=linewidth
            )
            ax.add_patch(ell)

        # 4) 绘图
        plt.figure(figsize=(4, 4), dpi=140)
        ax = plt.gca()

        palette = sns.color_palette(palette, n_colors=df["condition"].nunique())
        sns.scatterplot(
            data=df, x="PC1", y="PC2", hue="condition",
            s=30, edgecolor="white", linewidth=1.0, palette=palette, ax=ax
        )

        # 为每个分组叠加置信椭圆
        color_map = dict(zip(df["condition"].unique(), palette))
        for g in df["condition"].unique():
            sub = df[df["condition"] == g][["PC1", "PC2"]].values
            if len(sub) >= min_samples_for_ellipse and draw_ellipse:  # 检查条件
                plot_confidence_ellipse(sub, ax, edgecolor=color_map[g], facecolor=color_map[g], z=0.95)

        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()

        pad_x = 0.2 * (x_max - x_min)  # 20% 外边距
        pad_y = 0.2 * (y_max - y_min)

        ax.set_xlim(x_min - pad_x, x_max + pad_x)
        ax.set_ylim(y_min - pad_y, y_max + pad_y)

        ax.set_xlabel(f"PC1 ({expl[0] * 100:.1f}%)")
        ax.set_ylabel(f"PC2 ({expl[1] * 100:.1f}%)")
        ax.set_title("PCA with 95% Confidence Ellipses")

        plt.legend(title="Condition", frameon=False)
        plt.tight_layout()

        # 保存图像
        if save_to is not None:
            plt.savefig(save_to)
        
        plt.show()
    
    

    def plsda_plt(self, T_scores, y, save_to=None):
        """
        绘制 PLS 结果和决策区域并可选择保存图像。
        
        参数：
        - T_scores: PLS 变换后的得分 (n_samples, n_components) 的数组。
        - y: 分组标签 (n_samples,) 的数组或列表。
        - class_names: 分组的类别名称列表。
        - save_to: 保存图像的文件名，默认为 None（即不保存）。
        """

        from sklearn.linear_model import LogisticRegression

        class_names = np.unique(y)

        # 得分 DataFrame
        df_scores = pd.DataFrame(T_scores, columns=["LV1", "LV2"])
        df_scores["group"] = y

        # 1) 在得分空间上拟合一个简单分类器用于画决策边界
        clf = LogisticRegression(multi_class="auto", max_iter=1000)
        clf.fit(df_scores[["LV1", "LV2"]].to_numpy(), y)

        # 2) 创建网格
        pad = 0.10  # 边界留白比例
        x_min, x_max = df_scores["LV1"].min(), df_scores["LV1"].max()
        y_min, y_max = df_scores["LV2"].min(), df_scores["LV2"].max()
        x_pad = (x_max - x_min) * pad
        y_pad = (y_max - y_min) * pad
        x_min, x_max = x_min - x_pad, x_max + x_pad
        y_min, y_max = y_min - y_pad, y_max + y_pad

        xx, yy = np.meshgrid(
            np.linspace(x_min, x_max, 400),
            np.linspace(y_min, y_max, 400)
        )
        grid = np.c_[xx.ravel(), yy.ravel()]

        # 3) 预测网格类别（或概率）
        Z = clf.predict(grid)  # 类别标签
        Z = Z.reshape(xx.shape)

        # 4) 绘图
        plt.figure(figsize=(4, 4), dpi=140)
        ax = plt.gca()
        palette = sns.color_palette("Set2", n_colors=len(class_names))
        palette_map = {c: col for c, col in zip(class_names, palette)}

        # 决策区域底图（使用较浅的颜色）
        from matplotlib.colors import ListedColormap
        idx_map = {c: i for i, c in enumerate(class_names)}
        Z_idx = np.vectorize(idx_map.get)(Z)
        cmap_light = ListedColormap([sns.desaturate(palette_map[c], 0.9) for c in class_names])

        ax.contourf(xx, yy, Z_idx, alpha=0.25, cmap=cmap_light, levels=len(class_names))

        # 决策边界线（等概率/分界线）
        ax.contour(xx, yy, Z_idx, levels=np.arange(len(class_names)), colors="k", alpha=0.2, linewidths=0.8)

        # 叠加样本散点
        sns.scatterplot(
            data=df_scores, x="LV1", y="LV2", hue="group",
            s=30, edgecolor="none", linewidth=1.0,
            palette=palette, ax=ax
        )

        # 标注与样式
        ax.set_xlabel("LV1 (PLS-DA)")
        ax.set_ylabel("LV2 (PLS-DA)")
        ax.set_title("PLS-DA scores with decision regions")

        # 方框坐标轴
        for s in ax.spines.values():
            s.set_visible(True)
            s.set_linewidth(1.2)

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.legend(title="Group", frameon=False)
        plt.tight_layout()

        # 保存图像
        if save_to is not None:
            plt.savefig(save_to)
            print(f"The plot is saved to: {save_to}")
        
        plt.show()


    def swatch_colors(self, colors):
        # 设置色块的大小（1.5 cm）  
        block_size_cm = 2.5
        block_size_inch = block_size_cm / 2.54  # 转换为英寸（1 cm = 0.393701 in）  

        # 创建图形和轴  
        fig, ax = plt.subplots(figsize=(len(colors) * block_size_inch, block_size_inch))

        # 在每个色块中显示颜色  
        for i, color in enumerate(colors):  
            ax.add_patch(plt.Rectangle((i, 0), 1, 1, color=color))  

        # 设置轴的范围和标签  
        ax.set_xlim(0, len(colors))  
        ax.set_ylim(0, 1)  
        ax.set_xticks([i + 0.5 for i in range(len(colors))])  # 设置 x 轴的刻度位置  
        ax.set_xticklabels(colors)  # 设置 x 轴的标签为颜色代码  
        ax.set_yticks([])  # 隐藏 y 轴刻度  

        # 设置标题  
        plt.title('Color Swatches')  
        plt.show()

    def swatch_self_colors(self, colors):
        self.swatch_colors(self.colors[colors])

    def tic_rt_mz(self, df,
                  x,
                  y,
                  size,
                  color=None,
                  alpha=0.55,
                  shape=None,
                  palette='Nature_2_(5)'):  

        # 计算 y 轴的最小值和最大值  
        y_min = df[y].min() - 50  
        y_max = df[y].max() + 100  

        if palette in self.colors:
            palette = self.colors[palette]
        
        if color and shape:
            plot = (  
                ggplot(df, aes(x=x, y=y, size=size, color=color, shape=shape)) +  # 点大小依据 log_pkarea, 颜色依据 ionmode   
                geom_point(alpha=alpha)
                ) 
        elif color:
            plot = (  
                ggplot(df, aes(x=x, y=y, size=size, color=color)) +  # 点大小依据 log_pkarea, 颜色依据 ionmode   
                geom_point(alpha=alpha)
                )  
        elif shape:
            df[shape] = df[shape].astype(str)  
            plot = (  
                ggplot(df, aes(x=x, y=y, size=size, shape=shape)) +  # 点大小依据 log_pkarea, 颜色依据 ionmode   
                geom_point(color=self.colors[palette][0]) +
                scale_shape_manual(values={'False': 'o', 'True': '^'})  # 定义形状 
                )  
        else:
            plot = (  
                ggplot(df, aes(x=x, y=y, size=size)) +  # 点大小依据 log_pkarea, 颜色依据 ionmode   
                geom_point(color=self.colors[palette][0])
                )                                     

        # 绘制散点图  
        plot = (plot +             
            labs(x='Retention Time (min)', y='Precursor m/z', color='Ion Mode', size='Log (peak area)') +  
            ylim(y_min, y_max) +  # 设置纵坐标范围  
            scale_color_manual(values=palette) +
            theme_classic()+
            theme(aspect_ratio=2/3)  # 设置长宽比  
        )  

        return plot  


    def volcano(self, df, x, y, fill,
                xcut = 1, ycut = 2,
                title = '',
                xlab = r'$\mathrm{log_{2} \ Fold Change}$', # 使用LaTeX语法设置x轴标签 
                ylab = r'-$\log_{10}(\mathrm{p\text{-}value})$',# 使用LaTeX语法设置y轴标签
                palette='Set1',
                save_to=None):
        
        colors = plt.get_cmap(palette)([0, 1])
        colors = ['#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))\
                    for rgba in colors]

        x_limit = max(abs(min(df[x])), abs(max(df[x]))) 

        volcano_plot = (  
            ggplot(df, aes(x=x, y=y, color=fill)) +  
            geom_point (alpha=0.5, size=2, shape='o', stroke=0) +  # Set transparency for points  
            labs(title=f'{title}, Volcano Plot', x=xlab, y=ylab) +    
            xlim(-x_limit, x_limit) +
            ylim(0, max(2.5, df[y].max())) +
            scale_color_manual(values={
                                        'up': colors[0],
                                        'dn': colors[1],
                                         
                                        'no': '#D3D3D3' # light grey  
                                    }) +
            geom_vline(xintercept=-xcut, linetype='dashed', color='grey') +  
            geom_vline(xintercept= xcut, linetype='dashed', color='grey') + 
            geom_hline(yintercept= ycut, linetype='dashed', color='grey') +
            self.theme
        )  

        if save_to:
           ggsave(volcano_plot, save_to, dpi=300) 

        return volcano_plot

    def venn(self, data, palette='Set1', alpha=0.65, save_to:str = None):
        plot = Venn(data,
                    palette = palette,
                    fontsize = self.fontsize,
                    alpha=alpha,
                    save_to = save_to)
        return plot



class Venn:
    """
    A class to plot Venn diagrams for 2, 3, or 4 sets.
    Modified from pyvenn: https://github.com/tctianchi/pyvenn
    """

    def __init__(self, data, fill=['number'], palette='Set1', fontsize=14, alpha=0.65, save_to=None):
        """
        Initialize the Venn diagram plotter.

        Parameters:
            data (dict): A dictionary where keys are set names and values are lists of elements.
            fill (list): Options for labeling: ["number", "logic", "percent"].
            palette (str): Color palette name.
            fontsize (int): Font size for labels.
            alpha (float): Transparency of the circles.
            save_to (str): File path to save the plot. If None, plot is displayed.
        """
        plt.clf()  # Clear current figure
        self.fig = plt.figure(0, figsize=(9, 7), dpi=96)
        self.ax = self.fig.add_subplot(111, aspect='equal')
        self.ax.set_axis_off()
        self.ax.set_ylim(bottom=0.0, top=1.0)
        self.ax.set_xlim(left=0.0, right=1.0)

        self.palette = ListedColormap(plt.get_cmap(palette).colors)
        self.fontsize = fontsize
        self.alpha = alpha

        # Validate input data
        if not isinstance(data, dict):
            raise ValueError("Input data must be a dictionary.")
        if len(data) not in (2, 3, 4):
            raise ValueError("The data length must be 2, 3, or 4.")

        # Prepare data
        values = list(data.values())
        keys = list(data.keys())
        labels = self._divide(values, fill=fill)

        # Draw Venn diagram
        draw_func = getattr(self, f'_venn{len(data)}')
        draw_func(labels, keys)

        # Save or show the plot
        if save_to:
            self.fig.savefig(save_to, transparent=True)
        else:
            plt.show()

    def _divide(self, data, fill=["number"]):
        """
        Generate labels for the Venn diagram regions.

        Parameters:
            data (list): List of sets.
            fill (list): Options for labeling: ["number", "logic", "percent"].

        Returns:
            dict: A dictionary of labels for each region.
        """
        N = len(data)
        sets_data = [set(data[i]) for i in range(N)]  # Convert lists to sets
        s_all = set(chain(*data))  # Union of all sets
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

    def _draw_ellipse(self, xy, width, height, angle=0, color_index=0):
        """Draw an ellipse on the plot."""
        self.ax.add_patch(patches.Ellipse(
            xy=xy, width=width, height=height, angle=angle,
            alpha=self.alpha, color=self.palette(color_index)
        ))

    def _draw_text(self, x, y, text, ha='center', va='center'):
        """Draw text on the plot."""
        self.ax.text(x, y, text, fontsize=self.fontsize, ha=ha, va=va)

    def _venn2(self, labels, names):
        """Draw a 2-set Venn diagram."""
        self.ax.set_ylim(bottom=0.0, top=0.7)
        # Draw ellipses
        self._draw_ellipse((0.375, 0.3), 0.5, 0.5, color_index=0)
        self._draw_ellipse((0.625, 0.3), 0.5, 0.5, color_index=1)
        # Draw labels
        self._draw_text(0.74, 0.30, labels.get('01', ''))
        self._draw_text(0.26, 0.30, labels.get('10', ''))
        self._draw_text(0.50, 0.30, labels.get('11', ''))
        self._draw_text(0.20, 0.56, names[0])
        self._draw_text(0.80, 0.56, names[1])

    def _venn3(self, labels, names):
        """Draw a 3-set Venn diagram."""
        # Draw ellipses
        self._draw_ellipse((0.333, 0.633), 0.55, 0.55, color_index=0)
        self._draw_ellipse((0.666, 0.633), 0.55, 0.55, color_index=1)
        self._draw_ellipse((0.500, 0.310), 0.55, 0.55, color_index=2)
        # Draw labels
        self._draw_text(0.50, 0.27, labels.get('001', ''))
        self._draw_text(0.73, 0.65, labels.get('010', ''))
        self._draw_text(0.61, 0.46, labels.get('011', ''))
        self._draw_text(0.27, 0.65, labels.get('100', ''))
        self._draw_text(0.39, 0.46, labels.get('101', ''))
        self._draw_text(0.50, 0.65, labels.get('110', ''))
        self._draw_text(0.50, 0.51, labels.get('111', ''))
        self._draw_text(0.15, 0.87, names[0])
        self._draw_text(0.85, 0.87, names[1])
        self._draw_text(0.50, 0.02, names[2])

    def _venn4(self, labels, names):
        """Draw a 4-set Venn diagram."""
        o = 0.500  # Center of the plot
        dx = 0.18
        dy = 0.08
        # Draw ellipses
        self._draw_ellipse((o - dx, o - dy), 4 * dx, 2 * dx, angle=135, color_index=0)
        self._draw_ellipse((o, o), 4 * dx, 2 * dx, angle=135, color_index=1)
        self._draw_ellipse((o, o), 4 * dx, 2 * dx, angle=45, color_index=2)
        self._draw_ellipse((o + dx, o - dy), 4 * dx, 2 * dx, angle=45, color_index=3)
        # Draw labels
        label_positions = [
            (o + dx * 2.00, o + dy * 0.50, labels.get('0001', '')),
            (o + dx * 0.75, o + dy * 2.50, labels.get('0010', '')),
            (o + dx * 1.25, o + dy * 1.25, labels.get('0011', '')),
            (o - dx * 0.75, o + dy * 2.50, labels.get('0100', '')),
            (o + dx, o - dy * 2.00, labels.get('0101', '')),
            (o, o + dy * 1.25, labels.get('0110', '')),
            (o + dx * 0.75, o - dy * 0.25, labels.get('0111', '')),
            (o - dx * 2.00, o + dy * 0.50, labels.get('1000', '')),
            (o, o - dy * 3.75, labels.get('1001', '')),
            (o - dx, o - dy * 2.00, labels.get('1010', '')),
            (o - dx * 0.25, o - dy * 2.75, labels.get('1011', '')),
            (o - dx * 1.25, o + dy * 1.25, labels.get('1100', '')),
            (o + dx * 0.25, o - dy * 2.75, labels.get('1101', '')),
            (o - dx * 0.75, o - dy * 0.25, labels.get('1110', '')),
            (o, o - dy * 1.75, labels.get('1111', '')),
            (o - dx * 2.25, o + dy * 2.75, names[0]),
            (o - dx * 1.00, o + dy * 3.75, names[1]),
            (o + dx * 1.00, o + dy * 3.75, names[2]),
            (o + dx * 2.25, o + dy * 2.75, names[3])
        ]
        for x, y, text in label_positions:
            self._draw_text(x, y, text)
    
    def save(self, fname, dpi=300, transparent=True):
        return self.fig.savefig(fname=fname, dpi=dpi, transparent=transparent)
