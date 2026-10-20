# Frameshift P0 修复报告

- 日期：2026-10-20
- 仓库：D:\Documents\Codes\java\Frameshift
- 基线 commit：`98e8381`
- 最终 commit：`72fb455`
- 标签：`v0.1-marked` / `v0.2-marked` / `v0.3-p0-fixed`
- 工具：D:\Documents\Grenzen\jgrenzen.py

---

## 一、概览

本次工程分两阶段：

1. **块标记**：15 个 Java 文件全部加 `// GRENZEN` 块标记，共 71 个块（15 顶层类 + 4 内部类 + 8 内部类方法 + 44 方法），142 行标记。
2. **P0 修复**：修复"错误分支退出码恒为 0"这一影响所有自动化调用的缺陷，并解耦工具类的 `System.exit` 副作用。

修复规模：3 次提交，10 个文件，29 处代码变更。

---

## 二、体检方法

1. **编译层**：javac 24.0.1 全量编译 + `-Xlint:all`。
2. **静态模式扫描**：两轮，覆盖空 catch、String 字面量 `==`、死循环、资源泄漏、`System.exit` 位置、parse 无 try 等 13 类模式。
3. **人工定性**：30 处 `System.exit` 逐处看上下文，区分"错误分支"与"正常退出"。

---

## 三、编译层结论

- javac rc = 0，**0 错误**。
- 6 个警告，全部无害：
  - 3 个 `[cast]` 冗余 double 转换（Pearson / RandomCodes / SubScore）
  - 3 个 `[deprecation]` `Runtime.exec(String)`（Similarity / SubScore / shiftableExons3）

---

## 四、P0 缺陷

### 4.1 现象

30 处 `System.exit` 中，28 处位于错误分支（文件找不到、参数非法、格式错误、外部程序失败等），但**全部传 `0`**。

### 4.2 后果

- `System.exit(0)` 表示成功。错误分支返回 0，等于告诉操作系统"一切正常"。
- bash `&&` 链、批处理脚本、CI 流水线**无法感知失败**。
- 出错时产生的半截输出文件被当作正常产物。
- 返回码永远是 0，日志与监控看不见错误。

### 4.3 定性

**P0 功能缺陷**，不是风格问题。单次手动运行不受影响（人会看屏幕），但任何自动化调用都不可信。

---

## 五、修复明细

### 5.1 d1：错误分支 `exit(0)` -> `exit(1)`

commit `72fb455`，7 个文件，18 处替换。

| 文件 | 处数 | 触发条件 |
|---|---|---|
| FrameTranslation.main | 1 | args 数量 |
| MsGeneticCode.main | 3 | args 数量 / 参数范围 / 任务数超限 |
| RandomCDSs.main | 4 | 密码子数与序列数守卫 |
| RandomCodes.main | 1 | 参数组合 |
| Similarity.main | 3 | args / SeqType / clustalw 失败 |
| SubScore.main | 3 | 同上 |
| shiftableExons3 | 3 | usage / clustalw 失败 / blastp 头格式 |

### 5.2 d2：工具类抛 IOException

commit `abeedde`，3 个文件，9 处替换 + 1 处插入。

| 文件 | UUID | 方法 | 改动 |
|---|---|---|---|
| AAindex.java | `f4f0fb10` | readAAidx | 2x `exit(0)` -> `throw IOException` |
| AAindex.java | `73853b4d` | readAAidxTable | 3x 同上 |
| ScoringMatrix.java | `12ff9921` | readMatrix | 4x 同上 |
| ShiftCodons.java | `0a177bd0` | main | 1x `exit(0)`->`exit(1)`；catch 末尾插入 `exit(1)` |

三个工具方法原本已声明 `throws Exception`，改抛 `IOException` 无需改签名。

调用方处理：

- `MsGeneticCode.main` / `AlternativeCodes.main` / `RandomCodes.main` 声明 `throws Exception`，异常冒泡到 JVM，返回码 1，**零改动**。
- `ShiftCodons.main` 无 `throws`，但调用点在 try 内，catch 补 `exit(1)` 显式失败。

### 5.3 保留项

| 位置 | 原因 |
|---|---|
| Similarity.java:353 | 生成 MSA 命令文件后提示用户手动执行，属**正常退出**，保留 `exit(0)` |
| SubScore.java:353 | 同上 |

---

## 六、修复后验证

| 检查项 | 结果 |
|---|---|
| exit(0) 剩余 | 2 处（Similarity:353 / SubScore:353），均正常退出 |
| exit(1) 总数 | 20（18 替换 + 1 替换 + 1 插入） |
| 其他返回值 | 0 |
| javac 全量编译 | rc = 0 |
| 块结构完整性 | 71 块 BEGIN/END 配对，栈匹配，无嵌套错误 |

数字对账：`30 = 9 throw + 1 exit1(d2) + 18 exit1(d1) + 2 保留`。

---

## 七、git 历史

```
72fb455  d1: exit(0) -> exit(1) x18 / 7 files
abeedde  d2: IOException + catch exit(1)
bd986c6  Translation re-mark
2bb6d22  initial GRENZEN markers (71 blocks / 15 files)
98e8381  baseline
```

标签：

- `v0.1-marked`：首轮标记完成
- `v0.2-marked`：Translation 重标记 + jgrenzen 稳定
- `v0.3-p0-fixed`：P0 修复完成

任何时刻可精确回滚。

---

## 八、工具链

`D:\Documents\Grenzen\jgrenzen.py`，Java 版块标记工具，patchlist 格式与 grenzen.py 完全兼容。

```
jgrenzen.py scan <file> [--out PATH]
jgrenzen.py mark-blocks <file> [--dry-run] [--no-inner] [--no-methods]
jgrenzen.py check <patchlist>
jgrenzen.py apply <patchlist> <file> --dry-run
jgrenzen.py apply <patchlist> <file> --safe [--javac]
jgrenzen.py status <file>
```

- 标记格式：`// GRENZEN { <uuid> <name>` / `// } GRENZEN <uuid>`
- 支持内部类与内部类方法（三级结构）
- `--safe`：块配对自检 + 自动备份，不调 javac；`--javac` 单开
- 换行保持原文件风格，二进制写盘

---

## 九、残留演进债

| 优先级 | 项 | 规模 | 说明 |
|---|---|---|---|
| P1 | 巨型 main 拆分 | 4 个方法 2473 行 | Similarity.main 742 / SubScore.main 696 / MsGeneticCode.main 534 / ShiftCodons.main 501 |
| P2 | printStackTrace 日志化 | 7 处 | 风格问题，单机科研程序可接受 |
| P3 | 冗余 cast / deprecation 清理 | 6 处 | 无害 |
| P4 | 内部类命名小写开头 | 4 个类 | 改动需连带调用点，成本高收益低，建议不动 |

无 P0 存量。

---

## 十、后续建议

1. 当前代码可直接用于生产与论文复现。
2. 若长期维护，优先做 P1（巨型 main 拆分），但需单独立项、逐方法迁移。
3. 若扩展到 `AlternativeCodes/` 等其他目录的 Java 文件，直接复用 jgrenzen.py 的 mark-blocks 流程。
4. 任何后续改动走 `scan -> patchlist -> check -> dry-run -> safe` 五步，dry-run 全绿才写盘。
