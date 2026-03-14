# tmux_tut

## Quick Start

```bash
ssh -p 9211 <user>@220.248.74.126
tmux new-session -A -s main
```

## Core Commands

```bash
tmux new -s main
tmux ls
tmux attach -t main
tmux attach -d -t main
tmux new-session -A -s main
tmux kill-session -t main
tmux kill-server
tmux attach -d -t main || tmux new -s main
```

## Shortcuts (prefix: Ctrl-b)

```text
Detach: Ctrl-b d
New window: Ctrl-b c
Switch window: Ctrl-b n / Ctrl-b p / Ctrl-b 0..9
Split left-right: Ctrl-b %
Split up-down: Ctrl-b "
Switch pane: Ctrl-b + arrow
Close pane: exit
```

## Sessions

```bash
tmux new-session -A -s main
tmux new-session -A -s pipe
tmux new-session -A -s ref
tmux new-session -A -s monitor
```

## Snakemake Example

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/snakemake \
  -s /path/to/project/Snakefile \
  --configfile /path/to/project/build_ref.yaml \
  --cores 32 --rerun-incomplete -p reference_qc reference \
  2>&1 | tee /path/to/project/logs/ref_$(date +%F_%H%M%S).log
```

## Python Script Example

```bash
/biosoftware/miniconda/envs/snakemake_env/bin/python /path/to/project/scripts/bam_uniformity_qc.py \
  --target-bam /path/to/sample.bam \
  --ref-bams /path/to/ref1.bam /path/to/ref2.bam /path/to/ref3.bam \
  --bin-size 200000 \
  --mapq 30 \
  --outdir /path/to/project/qc/sample01 \
  --log /path/to/project/qc/sample01/bam_uniformity_qc.log
```

## Monitoring

```bash
tail -f /path/to/project/logs/wisecondorx/build_reference.log
htop
watch -n 2 "df -h; echo '---'; free -h; echo '---'; uptime"
```

## zsh aliases (~/.zshrc)

```bash
alias tmain='tmux new-session -A -s main'
alias tpipe='tmux new-session -A -s pipe'
alias tref='tmux new-session -A -s ref'
alias tmonitor='tmux new-session -A -s monitor'
```

```bash
source ~/.zshrc
```
