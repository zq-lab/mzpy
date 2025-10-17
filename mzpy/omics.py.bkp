'''
组学数据处理

-----------------------------------------
Chat: 
像代谢组或转录组的数据，至少是三个个有主键关联的数据表，一个数据表用于存储基因特征或MS特征信息，
另一个用于存储每个特征在每个样本中的定量结果，还有一个数据表需要存放每个样本编号的分组，
而进行两个组的对比后，还需要存储进行对比分析后的log2 fold change和校正p值.
请你写一个python的类定义，如何保存这样的多个表，并提供对比函数vs，vs接收两个组G1、G2，计算G1/G2的均值比值和p值
------------------------------------------

下面给出一个简洁实用的 Python 类，用于：

保存四类核心数据表：
特征表（基因/代谢物特征信息）
表达/定量矩阵（特征×样本）
样本分组表（样本→分组）
对比结果表（可按需生成与缓存）
提供 vs(G1, G2) 方法，计算两组间：
每个特征的均值比值（以 log2 fold change 表达：log2(mean(G1)/mean(G2))）
显著性检验 p 值（并给出 Benjamini-Hochberg FDR 校正后的 q 值）

'''

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Literal, Dict
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


@dataclass
class ContrastResult:
    """
    存放单次对比的结果与元信息。
    """
    lfc: pd.Series               # log2 fold change: log2(mean(numerator)/mean(denominator))
    pval: pd.Series              # 未校正 p 值
    qval: pd.Series              # 多重校正后的 p 值（FDR）
    mean_numerator: pd.Series    # 分子组的均值（线性尺度）
    mean_denominator: pd.Series  # 分母组的均值（线性尺度）
    n_numerator: int             # 分子组样本数
    n_denominator: int           # 分母组样本数
    method: str                  # 统计检验方法说明
    name: str                    # 对比名称，如 "Treatment_vs_Control"
    features: Optional[pd.DataFrame] = None  # 可选：保留特征注释子集

    def to_dataframe(self, include_features: bool = True) -> pd.DataFrame:
        df = pd.DataFrame({
            "mean_numerator": self.mean_numerator,
            "mean_denominator": self.mean_denominator,
            "log2FC": self.lfc,
            "pval": self.pval,
            "qval": self.qval,
        })
        if include_features and self.features is not None:
            df = df.join(self.features, how="left")
        return df


class OmicsDataset:
    """
    管理组学数据的多表结构，并进行两组对比分析。

    数据对象：
    - features: 行为特征（feature_id 唯一），列为特征注释信息
    - matrix: 行为特征，列为样本 ID，值为定量结果（线性尺度）
    - sample_groups: 索引为样本 ID，至少包含一列 'group' 指明分组
    """

    def __init__(
        self,
        features: pd.DataFrame,
        matrix: pd.DataFrame,
        sample_groups: pd.DataFrame,
        feature_id_col: str = "feature_id",
        group_col: str = "group",
    ):
        # 1) 处理 features：保证 feature_id 为索引
        if feature_id_col not in features.columns and features.index.name != feature_id_col:
            raise ValueError(f"features 缺少主键列：{feature_id_col}")
        if features.index.name != feature_id_col:
            features = features.set_index(feature_id_col)
        if features.index.duplicated().any():
            dup = features.index[features.index.duplicated()].unique().tolist()[:5]
            raise ValueError(f"features 存在重复 feature_id，如: {dup} ...")

        # 2) 处理 matrix：行映射到 feature_id，列映射到样本
        if matrix.index.name != features.index.name and features.index.name in matrix.columns:
            matrix = matrix.set_index(features.index.name)

        # 确保数值类型
        matrix = matrix.apply(pd.to_numeric, errors="coerce")

        # 3) 处理 sample_groups：保证样本 ID 为索引，且包含 group 列
        if sample_groups.index.name is None and "sample_id" in sample_groups.columns:
            sample_groups = sample_groups.set_index("sample_id")
        if group_col not in sample_groups.columns:
            raise ValueError(f"sample_groups 缺少分组列：{group_col}")
        if sample_groups.index.duplicated().any():
            dup = sample_groups.index[sample_groups.index.duplicated()].unique().tolist()[:5]
            raise ValueError(f"sample_groups 存在重复样本 ID，如: {dup} ...")

        # 4) 对齐样本：matrix 列与 sample_groups 索引取交集
        common_samples = matrix.columns.intersection(sample_groups.index)
        if len(common_samples) == 0:
            raise ValueError("matrix 与 sample_groups 无共同样本，无法对齐。")
        if len(common_samples) < matrix.shape[1]:
            matrix = matrix.loc[:, common_samples]
        sample_groups = sample_groups.loc[common_samples]

        # 5) 对齐特征：features 与 matrix 行取交集
        common_features = matrix.index.intersection(features.index)
        if len(common_features) == 0:
            raise ValueError("features 与 matrix 无共同特征，无法对齐。")
        if len(common_features) < matrix.shape[0]:
            matrix = matrix.loc[common_features]
        features = features.loc[common_features]

        # 存储
        self.features = features
        self.matrix = matrix
        self.sample_groups = sample_groups
        self.feature_id_col = feature_id_col
        self.group_col = group_col

        self.contrasts: Dict[str, ContrastResult] = {}

    # ---------------------------
    # 缺失/零值处理
    # ---------------------------
    def impute_min_fraction(
        self,
        *,
        inplace: bool = True,
        return_value: bool = False,
    ) -> Optional[float]:
        """
        将矩阵中的 0 值与 NaN 替换为：全矩阵正数非零最小值的 1/10。
        仅当存在需要替换的值（NaN、<=0）时才执行替换；否则不做任何处理。

        规则：
        - 仅视 >0 的值为“正数非零”。<=0（如 0 或负数）与 NaN 视为待替换。
        - 若全矩阵不存在正数非零值且存在待替换值，则抛出 ValueError（避免引入武断尺度）。

        参数：
        - inplace: True 则原地修改 self.matrix；False 则不修改，仅返回替换值（若 return_value=True）。
        - return_value: 是否返回用于替换的数值（便于记录）。仅在计算出替换值时返回。

        返回：
        - 若 return_value=True 且执行/可执行替换，则返回替换用的数值（float）；否则返回 None。
        """
        # 将 DataFrame 视作数组便于快速检测
        arr = self.matrix.to_numpy(dtype=float, copy=False)

        # 检测是否存在需要替换的位置（NaN 或 <= 0）
        mask_nan = np.isnan(arr)
        mask_nonpos = arr <= 0  # 包含 0 与 负数（若有）
        need_replace = np.any(mask_nan | mask_nonpos)

        if not need_replace:
            # 不存在 0 或 NaN（或负值），直接返回
            if return_value:
                return None
            return None

        # 找到正数非零的最小值
        mask_valid = ~np.isnan(arr)
        mask_pos = (arr > 0) & mask_valid
        if not np.any(mask_pos):
            # 有需要替换的值，但没有任何正数可作为参照
            raise ValueError("矩阵中没有正数非零值，无法根据规则确定替换值。")

        nonzero_min = float(np.min(arr[mask_pos]))
        replace_value = nonzero_min / 10.0

        if inplace:
            # 使用 DataFrame 的 where 进行矢量化替换：
            # where(condition, other=...)：保留 condition 为 True 的原值，否则替换为 other
            mask_replace_df = self.matrix.isna() | (self.matrix <= 0)
            self.matrix = self.matrix.where(~mask_replace_df, other=replace_value)

        if return_value:
            return replace_value
        return None

    # ---------------------------
    # 实用工具
    # ---------------------------
    def _get_group_samples(self, group_name: str) -> pd.Index:
        samples = self.sample_groups.index[self.sample_groups[self.group_col] == group_name]
        if len(samples) == 0:
            raise ValueError(f"分组 {group_name} 中没有样本。")
        return samples

    # ---------------------------
    # 统计对比：vs(numerator_group, denominator_group)
    # ---------------------------
    def vs(
        self,
        numerator_group: str,
        denominator_group: str,
        *,
        test: Literal["welch", "ttest_ind"] = "welch",
        equal_var: Optional[bool] = None,
        tail: Literal["two-sided", "greater", "less"] = "two-sided",
        fdr_method: str = "fdr_bh",
        min_non_na: int = 2,
        cache: bool = True,
        keep_features_in_result: bool = True,
    ) -> ContrastResult:
        """
        计算 log2FC = log2(mean(numerator)/mean(denominator))，
        以及每个特征的 p 值与 FDR 校正 q 值，并返回 ContrastResult。

        注意：
        - 本函数不进行零值/空值替换。
        - 若你的数据包含 0 或 NaN，请在调用前使用 impute_min_fraction() 进行替换。
        """
        # 样本集合
        s_num = self._get_group_samples(numerator_group)
        s_den = self._get_group_samples(denominator_group)

        # 组别子矩阵（按列）
        X_num = self.matrix.loc[:, self.matrix.columns.intersection(s_num)]
        X_den = self.matrix.loc[:, self.matrix.columns.intersection(s_den)]

        # 计算均值（线性尺度）
        mean_num = X_num.mean(axis=1, skipna=True)
        mean_den = X_den.mean(axis=1, skipna=True)

        # log2FC（若 mean_den 为 0 或 NaN，可能出现 inf/NaN）
        with np.errstate(divide="ignore", invalid="ignore"):
            lfc = np.log2(mean_num / mean_den)

        # 选择 t 检验函数
        if test == "welch":
            ttest = lambda a, b: stats.ttest_ind(a, b, equal_var=False, nan_policy="omit", alternative=tail)
        elif test == "ttest_ind":
            if equal_var is None:
                equal_var = True
            ttest = lambda a, b: stats.ttest_ind(a, b, equal_var=equal_var, nan_policy="omit", alternative=tail)
        else:
            raise ValueError("test 仅支持 'welch' 或 'ttest_ind'")

        # 逐特征检验
        pvals = pd.Series(index=self.matrix.index, dtype=float)
        n_num = X_num.shape[1]
        n_den = X_den.shape[1]
        for fid in self.matrix.index:
            a = X_num.loc[fid].to_numpy()
            b = X_den.loc[fid].to_numpy()
            # 每组最少非 NA 数量限制
            if np.sum(~np.isnan(a)) < min_non_na or np.sum(~np.isnan(b)) < min_non_na:
                pvals.loc[fid] = np.nan
                continue
            res = ttest(a, b)
            pvals.loc[fid] = float(res.pvalue)

        # 多重校正（BH-FDR）
        valid = pvals.notna()
        qvals = pd.Series(index=pvals.index, dtype=float)
        if valid.any():
            _, q, _, _ = multipletests(pvals[valid].values, method=fdr_method)
            qvals.loc[valid] = q
        else:
            qvals.loc[:] = np.nan

        contrast_name = f"{numerator_group}_vs_{denominator_group}"

        result = ContrastResult(
            lfc=lfc,
            pval=pvals,
            qval=qvals,
            mean_numerator=mean_num,
            mean_denominator=mean_den,
            n_numerator=n_num,
            n_denominator=n_den,
            method=f"{test} t-test ({tail}), FDR: {fdr_method}",
            name=contrast_name,
            features=self.features.copy() if keep_features_in_result else None,
        )

        if cache:
            self.contrasts[contrast_name] = result

        return result

    # ---------------------------
    # 便捷方法
    # ---------------------------
    def get_contrast(self, name: str) -> Optional[ContrastResult]:
        return self.contrasts.get(name)

    def export_contrast(self, name: str, path: str, include_features: bool = True) -> None:
        res = self.get_contrast(name)
        if res is None:
            raise KeyError(f"未找到对比结果：{name}")
        df = res.to_dataframe(include_features=include_features)
        df.to_csv(path, index=True)

    def groups(self) -> pd.Series:
        return self.sample_groups[self.group_col]

    def description(self) -> str:
        return (
            f"Features: {self.features.shape[0]} x {self.features.shape[1]} "
            f"| Matrix: {self.matrix.shape[0]} x {self.matrix.shape[1]} "
            f"| Groups: {self.sample_groups[self.group_col].nunique()} levels, "
            f"{self.sample_groups.shape[0]} samples"
        )

