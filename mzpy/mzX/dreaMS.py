import sys
from typing import List, Tuple

# ===== 环境检查 =====
def check_required_packages() -> Tuple[bool, List[str]]:
    """
    检查是否安装了所需的包
    
    Returns:
        (all_ok, missing_packages): 
            all_ok: 布尔值，True 表示所有包都已安装
            missing_packages: 缺失的包列表
    """
    required_packages = {
        'torch': 'PyTorch (torch)',
        'dreams': 'DreaMS (dreams)',
        'numpy': 'NumPy (numpy)',
        'pandas': 'Pandas (pandas)',
        'tqdm': 'TQDM (tqdm)',
        'h5py': 'HDF5 for Python (h5py)'
    }
    
    missing = []
    
    for package_name, display_name in required_packages.items():
        try:
            __import__(package_name)
        except ImportError:
            missing.append((package_name, display_name))
    
    return len(missing) == 0, missing


def print_installation_guide(missing_packages: List[Tuple[str, str]]):
    """打印安装指南"""
    
    print("\n" + "=" * 70)
    print("❌ 缺失必要的依赖包")
    print("=" * 70)
    
    print("\n缺失的包：")
    for package_name, display_name in missing_packages:
        print(f"  • {display_name} ({package_name})")
    
    print("\n📦 安装方法：\n")
    
    # PyTorch 特殊处理
    torch_missing = any(pkg[0] == 'torch' for pkg in missing_packages)
    if torch_missing:
        print("1️⃣ 安装 PyTorch (选择适合你的版本):")
        print("   • GPU 版本 (CUDA 12.1):")
        print("     pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121")
        print("   • GPU 版本 (CUDA 11.8):")
        print("     pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118")
        print("   • CPU 版本:")
        print("     pip install torch torchvision torchaudio")
        print()
    
    # 其他包
    other_packages = [pkg for pkg in missing_packages if pkg[0] != 'torch']
    if other_packages:
        print("2️⃣ 安装其他必需包:")
        pip_packages = ' '.join([pkg[0] for pkg in other_packages])
        print(f"   pip install {pip_packages}")
        print()
    
    # Dreams 特殊说明
    dreams_missing = any(pkg[0] == 'dreams' for pkg in missing_packages)
    if dreams_missing:
        print("3️⃣ 安装 DreaMS:")
        print("   pip install dreams")
        print("   或从源码安装:")
        print("   git clone https://github.com/yourusername/dreams.git")
        print("   cd dreams")
        print("   pip install -e .")
        print()
    
    print("=" * 70 + "\n")


# 在模块导入前进行检查
all_ok, missing = check_required_packages()

if not all_ok:
    print_installation_guide(missing)
    sys.exit(1)

# ==========================================================================

import numpy as np
import torch
import time
from typing import List, Union
from pathlib import Path
from tqdm import tqdm
from dreams.api import PreTrainedModel
from dreams.definitions import DREAMS_EMBEDDING
import dreams.utils.data as du
import dreams.utils.dformats as dformats
from torch.utils.data import DataLoader, Dataset


class RawSpectraDatasetMemory(Dataset):
    """直接在内存中处理原始质谱数据的Dataset，无需临时文件"""
    
    def __init__(self, spectra, prec_mzs, spec_preproc: du.SpectrumPreprocessor):
        """
        Args:
            spectra: 质谱列表，每个元素是 (N_peaks, 2) 的 numpy 数组
            prec_mzs: 前体 m/z 数组
            spec_preproc: 预处理器
        """
        self.spectra = spectra
        self.prec_mzs = prec_mzs
        self.spec_preproc = spec_preproc
    
    def __len__(self):
        return len(self.spectra)
    
    def __getitem__(self, i):
        """
        获取单个样本
        
        Args:
            i: 样本索引
        
        Returns:
            字典，包含预处理后的质谱和前体 m/z
        """
        spectrum = self.spec_preproc(
            self.spectra[i], 
            prec_mz=self.prec_mzs[i], 
            high_form=True  # ✅ 关键：输入格式是 (N_peaks, 2)
        )
        return {
            du.SPECTRUM: spectrum, 
            du.PRECURSOR_MZ: self.prec_mzs[i]
        }


class DreamsEmbedder:
    """
    DreaMS 嵌入生成器（无临时文件版本）
    
    特点：
    - 直接在内存中处理数据
    - 无需临时文件
    - 模型缓存，后续调用高效
    - 支持 GPU 和 CPU
    - 自动批处理
    
    使用示例：
        embedder = DreamsEmbedder(device='cuda')
        embedder.load_model()
        embs = embedder.generate_embeddings(prec_mz, msms)
        embedder.unload_model()
    """
    
    def __init__(self, device: Union[str, torch.device] = None):
        """
        初始化 DreamsEmbedder
        
        Args:
            device: 计算设备，可选 'cuda'、'cpu'、torch.device 对象
                   如果为 None，自动选择（优先 GPU）
        """
        if device is None:
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        else:
            self.device = torch.device(device) if isinstance(device, str) else device
        
        self.model = None
        self.loaded = False
        print(f"DreamsEmbedder initialized with device: {self.device}")
    
    def load_model(self, verbose: bool = True):
        """
        加载预训练的 DreaMS 模型
        
        只有第一次调用才会真正加载模型，后续调用直接返回
        
        Args:
            verbose: 是否打印加载过程
        """
        if self.loaded:
            if verbose:
                print("✓ Model already loaded")
            return
        
        if verbose:
            print(f"Loading DreaMS pre-trained model...")
        
        start_time = time.time()
        
        try:
            # 加载预训练模型
            pretrained_model = PreTrainedModel.from_name(DREAMS_EMBEDDING)
            self.model = pretrained_model.model
            
            # 转移到指定设备
            self.model = self.model.to(self.device)
            self.model.eval()  # 设置为评估模式
            
            load_time = time.time() - start_time
            self.loaded = True
            
            if verbose:
                print(f"✓ Model loaded successfully in {load_time:.2f}s")
                if self.device.type == 'cuda':
                    print(f"  GPU: {torch.cuda.get_device_name(0)}")
                    print(f"  VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
        
        except Exception as e:
            print(f"✗ Failed to load model: {e}")
            raise
    
    @staticmethod
    def _prepare_spectra(
        precursor_mz: Union[np.ndarray, List[float]],
        msms: Union[List[List[List[float]]], List[np.ndarray]]
    ) -> tuple:
        """
        准备质谱数据（转换为标准格式）
        
        Args:
            precursor_mz: 前体 m/z 数组
            msms: 二碎片质谱集合
        
        Returns:
            (precursor_mz_array, msms_list) 元组
        
        Raises:
            ValueError: 如果输入数据格式不正确
        """
        precursor_mz = np.asarray(precursor_mz, dtype=np.float32)
        
        if not isinstance(msms, list):
            raise ValueError("msms must be a list of spectra")
        
        msms_converted = []
        
        for i, spec in enumerate(msms):
            # 转换为 numpy 数组
            if isinstance(spec, np.ndarray):
                spec = spec.astype(np.float32)
            elif isinstance(spec, list):
                spec = np.array(spec, dtype=np.float32)
            else:
                raise ValueError(f"Unsupported spectrum format at index {i}: {type(spec)}")
            
            # 确保是 2D 数组
            if spec.ndim != 2:
                raise ValueError(
                    f"Spectrum {i} must be 2D, got {spec.ndim}D with shape {spec.shape}"
                )
            
            # 检查格式：如果是 (2, N_peaks)，转置为 (N_peaks, 2)
            if spec.shape[0] == 2 and spec.shape[1] != 2:
                spec = spec.T
            
            # 验证最终格式
            if spec.shape[1] != 2:
                raise ValueError(
                    f"Spectrum {i} must have shape (N_peaks, 2), got {spec.shape}"
                )
            
            msms_converted.append(spec)
        
        return precursor_mz, msms_converted
    
    def generate_embeddings(
        self,
        precursor_mz: Union[np.ndarray, List[float]],
        msms: Union[List[List[List[float]]], List[np.ndarray]],
        batch_size: int = 32,
        progress_bar: bool = True,
        n_highest_peaks: int = 100,
        return_numpy: bool = True,
        verbose: bool = True
    ) -> Union[np.ndarray, torch.Tensor]:
        """
        生成质谱数据的 DreaMS 嵌入向量
        
        Args:
            precursor_mz: 前体 m/z 数组或列表，形状 (n_spectra,)
            msms: 二碎片质谱集合，格式: [[[m/z1, intens1], [m/z2, intens2], ...], ...]
            batch_size: 批处理大小，默认 32
            progress_bar: 是否显示进度条，默认 True
            n_highest_peaks: 使用最高强度的峰数，默认 100
            return_numpy: 是否返回 numpy 数组，默认 True
            verbose: 是否打印处理过程，默认 True
        
        Returns:
            embeddings: numpy 数组或 torch tensor，形状 (n_spectra, 1024)
        
        Raises:
            RuntimeError: 如果模型未加载
            ValueError: 如果输入数据格式不正确
        """
        if not self.loaded:
            raise RuntimeError("Model not loaded. Call load_model() first.")
        
        start_time = time.time()
        
        # ===== 步骤 1：准备数据 =====
        if verbose:
            print("Preparing spectra...")
        
        precursor_mz, msms_converted = self._prepare_spectra(precursor_mz, msms)
        n_spectra = len(msms_converted)
        
        if verbose:
            print(f"✓ Prepared {n_spectra} spectra")
            print(f"  First spectrum shape: {msms_converted[0].shape}")
        
        # ===== 步骤 2：初始化预处理器 =====
        if verbose:
            print("Initializing preprocessor...")
        
        spec_preproc = du.SpectrumPreprocessor(
            dformat=dformats.DataFormatA(),
            n_highest_peaks=n_highest_peaks
        )
        
        # ===== 步骤 3：创建内存 Dataset =====
        if verbose:
            print("Creating dataset...")
        
        dataset = RawSpectraDatasetMemory(
            spectra=msms_converted,
            prec_mzs=precursor_mz,
            spec_preproc=spec_preproc
        )
        
        dataloader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=False,
            drop_last=False,
            num_workers=0  # 设为 0 避免多进程问题
        )
        
        # ===== 步骤 4：批量推理 =====
        if verbose:
            print(f"Generating embeddings for {n_spectra} spectra...")
        
        all_embeddings = []
        
        with torch.inference_mode():
            pbar = tqdm(
                dataloader,
                total=len(dataset),
                desc='Generating embeddings',
                disable=not progress_bar,
                unit='spectra'
            ) if progress_bar else dataloader
            
            for batch in pbar:
                batch_spectra = batch[du.SPECTRUM].to(self.device)
                
                # 推理
                embeddings = self.model(batch_spectra)
                
                # 转移到 CPU 并转换为 numpy
                embeddings = embeddings.cpu().detach().numpy()
                all_embeddings.append(embeddings)
        
        # ===== 步骤 5：合并结果 =====
        embeddings_array = np.concatenate(all_embeddings, axis=0)
        
        # 统计信息
        elapsed_time = time.time() - start_time
        throughput = n_spectra / elapsed_time
        
        if verbose:
            print(f"✓ Generated embeddings shape: {embeddings_array.shape}")
            print(f"  Total time: {elapsed_time:.2f}s")
            print(f"  Throughput: {throughput:.1f} spectra/sec")
        
        if return_numpy:
            return embeddings_array
        else:
            return torch.from_numpy(embeddings_array)
    
    def unload_model(self, verbose: bool = True):
        """
        卸载模型并释放 GPU 显存
        
        Args:
            verbose: 是否打印卸载过程
        """
        if self.model is not None:
            del self.model
            self.model = None
            self.loaded = False
            
            # 清空 GPU 缓存
            if self.device.type == 'cuda':
                torch.cuda.empty_cache()
            
            if verbose:
                print("✓ Model unloaded and memory freed")
        else:
            if verbose:
                print("✓ Model not loaded, nothing to unload")
    
    def get_model_info(self) -> dict:
        """获取模型信息"""
        info = {
            'loaded': self.loaded,
            'device': str(self.device),
            'model': type(self.model).__name__ if self.model else None,
        }
        
        if self.device.type == 'cuda':
            info['gpu_name'] = torch.cuda.get_device_name(0)
            info['gpu_memory_gb'] = torch.cuda.get_device_properties(0).total_memory / 1e9
        
        return info
    
    def __repr__(self) -> str:
        """字符串表示"""
        status = "loaded" if self.loaded else "unloaded"
        return f"DreamsEmbedder(device={self.device}, status={status})"

    def unload_model(self, verbose: bool = True):
        """旧接口：尽可能清理模型，但保持对外兼容性。
        说明：这方法不一定能完全释放显存，推荐使用 close/close_model 来彻底卸载。
        """
        try:
            if self.model is None:
                if verbose:
                    print("✔ No model to unload.")
                self.loaded = False
                return
            # 尝试把模型移出设备，避免显存锁定
            if hasattr(self.model, "cpu"):
                self.model = self.model.cpu()
            del self.model
            self.model = None
            self.loaded = False
            if verbose:
                print("✔ Model unloaded (best-effort).")
            torch.cuda.empty_cache()
        except Exception as e:
            if verbose:
                print(f"✗ Failed to unload model: {e}")

    def close(self, verbose: bool = True):
        """彻底关闭 embedder，卸载模型并释放 GPU 显存。
        等价于“卸载 + 释放显存 + 断开引用”。
        """
        if verbose:
            print("Closing DreamsEmbedder: unloading model and freeing resources...")

        # 1) 先尝试卸载模型
        if self.model is not None:
            try:
                # 将模型转到 CPU（如果在 CUDA 上）
                if getattr(self.model, "to", None) is not None:
                    try:
                        self.model = self.model.to('cpu')
                    except Exception:
                        pass
                # 删除引用
                del self.model
            except Exception as e:
                if verbose:
                    print(f"✗ Error while unloading model: {e}")

        self.model = None
        self.loaded = False

        # 2) 释放显存
        try:
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
                if verbose:
                    try:
                        if torch.cuda.is_available():
                            total = torch.cuda.get_device_properties(0).total_memory / 1e9
                            free = torch.cuda.memory_reserved(0) / 1e9
                            print(f"GPU memory cache cleared. (Total GPU memory: {total:.1f} GB, Reserved: {free:.2f} GB)")
                    except Exception:
                        pass
        except Exception as e:
            if verbose:
                print(f"✗ Failed to empty CUDA cache: {e}")

        if verbose:
            print("DreamsEmbedder closed.")


# ===== 使用示例 =====
if __name__ == '__main__':
    from peak import read_mgf
    import numpy as np
    
    # ===== 初始化嵌入器 =====
    print("=" * 70)
    print("Initialize DreamsEmbedder (No Temporary Files)")
    print("=" * 70)
    embedder = DreamsEmbedder(device='cuda')
    print(embedder)
    print()
    
    # ===== 加载模型 =====
    print("=" * 70)
    print("Load model")
    print("=" * 70)
    embedder.load_model()
    print("\nModel info:")
    for key, value in embedder.get_model_info().items():
        print(f"  {key}: {value}")
    print()
    
    # ===== 读取数据 =====
    df = read_mgf('../data/example_5_spectra.mgf')
    data = [ms for ms in df['MSMS'].values]
    
    # ===== 第一次生成嵌入向量 =====
    print("=" * 70)
    print("First call: Generate embeddings")
    print("=" * 70)
    import time
    start = time.time()
    embs1 = embedder.generate_embeddings(
        df['precursormz'].values, 
        data,
        verbose=True
    )
    time1 = time.time() - start
    print(f"Result shape: {embs1.shape}")
    print(f"First 5 values: {embs1[0, 300:305]}\n")
    
    # ===== 第二次生成嵌入向量（重用模型） =====
    print("=" * 70)
    print("Second call: Generate embeddings (using cached model)")
    print("=" * 70)
    start = time.time()
    embs2 = embedder.generate_embeddings(
        df['precursormz'].values, 
        data,
        verbose=True
    )
    time2 = time.time() - start
    print(f"Result shape: {embs2.shape}")
    print(f"First 5 values: {embs2[0, 300:305]}\n")
    
    # ===== 第三次生成嵌入向量（重用模型） =====
    print("=" * 70)
    print("Third call: Generate embeddings (using cached model)")
    print("=" * 70)
    start = time.time()
    embs3 = embedder.generate_embeddings(
        df['precursormz'].values, 
        data,
        verbose=False,
        progress_bar=True
    )
    time3 = time.time() - start
    print(f"Result shape: {embs3.shape}\n")
    
    # ===== 性能对比 =====
    print("=" * 70)
    print("Performance Comparison")
    print("=" * 70)
    print(f"First call:  {time1:.2f}s")
    print(f"Second call: {time2:.2f}s")
    print(f"Third call:  {time3:.2f}s")
    print(f"Speedup: {time1 / time2:.1f}x faster after caching\n")
    
    # ===== 验证结果一致性 =====
    print("=" * 70)
    print("Verify Consistency")
    print("=" * 70)
    print(f"embs1 == embs2? {np.allclose(embs1, embs2)}")
    print(f"embs2 == embs3? {np.allclose(embs2, embs3)}")
    print(f"Max difference: {np.max(np.abs(embs1 - embs2)):.2e}\n")
    
    # ===== 卸载模型 =====
    print("=" * 70)
    print("Unload model")
    print("=" * 70)
    embedder.unload_model()
    print()


### 相似度算法