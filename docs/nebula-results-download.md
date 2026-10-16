# Download and Synchronize Nebula Results

This is the canonical procedure for downloading data and results produced by
public analytical workflow jobs on Nebula. Run these commands in a terminal on the local
computer, not inside the Nebula SSH session. Connect to the UB VPN first.

## Canonical AI response

When asked how to download or synchronize HPC-produced project data, use this
procedure. Read the ignored `hpc/hpc_access.local.md` file for the SSH alias and
current HPC checkout path. Do not print credentials, private keys, passphrases,
IP addresses, or other secret connection values. If the local note is absent,
use `<NEBULA_SSH_ALIAS>` and `<HPC_REPO>` placeholders and ask the user to verify
them.

The examples below assume the SSH alias is `nebula`. Set `HPC_REPO` to the
checkout path recorded in the ignored local access notes.

## 1. Verify the remote checkout

From the local computer:

```bash
ssh nebula 'cd <HPC_REPO> && pwd && du -sh data/processed results logs 2>/dev/null'
```

Replace `<HPC_REPO>` before running the command. The output must show the
expected public analytical checkout. Stop if the path is wrong.

## 2. Enter the local checkout

```bash
cd /path/to/local/sagarmatha-microbiome-analysis
HPC_REPO=/path/to/nebula/sagarmatha-microbiome-analysis
```

## 3. Preview the synchronization

These dry runs do not change local files:

```bash
rsync -avhn --itemize-changes \
  nebula:"$HPC_REPO/data/processed/" \
  data/processed/

rsync -avhn --itemize-changes \
  nebula:"$HPC_REPO/results/" \
  results/

rsync -avhn --itemize-changes \
  nebula:"$HPC_REPO/logs/" \
  logs/
```

Review all proposed changes before continuing. In particular, confirm that the
source is Nebula and the destination is the local checkout.

For a Random Forest subject-validation run, prefer synchronizing only that
uniquely named run and its job log so unrelated local results are not touched:

```bash
RUN_ID=production_<PBS_JOB_ID>

rsync -avhn --itemize-changes \
  nebula:"$HPC_REPO/results/random_forest_subject_validation/$RUN_ID/" \
  "results/random_forest_subject_validation/$RUN_ID/"

rsync -avhn --itemize-changes \
  nebula:"$HPC_REPO/logs/random_forest_subject_validation.$RUN_ID.log" \
  logs/
```

## 4. Download with replacement backups

The following commands do not use `--delete`. Local-only files remain in place.
When an incoming file replaces a local file, rsync saves the previous local
version in a timestamped directory beside the repository.

```bash
STAMP=$(date +%Y%m%d_%H%M%S)
BACKUP_DIR="../sagarmatha-microbiome-analysis-before-nebula-sync-$STAMP"

mkdir -p "$BACKUP_DIR"

rsync -avh --partial --info=progress2 \
  --backup --backup-dir="$BACKUP_DIR/data-processed" \
  nebula:"$HPC_REPO/data/processed/" \
  data/processed/

rsync -avh --partial --info=progress2 \
  --backup --backup-dir="$BACKUP_DIR/results" \
  nebula:"$HPC_REPO/results/" \
  results/

rsync -avh --partial --info=progress2 \
  --backup --backup-dir="$BACKUP_DIR/logs" \
  nebula:"$HPC_REPO/logs/" \
  logs/
```

This synchronizes:

- generated datasets under `data/processed/`;
- all analysis outputs under `results/`;
- PBS and workflow logs under `logs/`.

The disposable `tmp/` directory and raw input data under `data/raw/` are not
downloaded by this procedure.

For the targeted Random Forest run, use a run-specific backup and keep results
and logs separate:

```bash
STAMP=$(date +%Y%m%d_%H%M%S)
RUN_ID=production_<PBS_JOB_ID>
BACKUP_DIR="../sagarmatha-microbiome-analysis-before-nebula-sync-$STAMP"
mkdir -p "$BACKUP_DIR/random-forest-subject-validation" "$BACKUP_DIR/logs"

rsync -avh --partial --info=progress2 \
  --backup --backup-dir="$BACKUP_DIR/random-forest-subject-validation" \
  nebula:"$HPC_REPO/results/random_forest_subject_validation/$RUN_ID/" \
  "results/random_forest_subject_validation/$RUN_ID/"

rsync -avh --partial --info=progress2 \
  --backup --backup-dir="$BACKUP_DIR/logs" \
  nebula:"$HPC_REPO/logs/random_forest_subject_validation.$RUN_ID.log" \
  logs/
```

Neither the general nor targeted commands use `--delete`.

## 5. Verify the local copy

```bash
du -sh data/processed results logs
git status --short --ignored
```

Generated data, results, logs, and archives are ignored by Git intentionally.
Do not force-add them with `git add -f`. Preserve final artifacts using the
local copy and an appropriate independent backup or publication archive.

## Troubleshooting

- `Could not resolve hostname nebula`: confirm the `nebula` entry exists in
  `~/.ssh/config` and that the UB VPN is connected.
- `No such file or directory`: verify `HPC_REPO` using the remote-checkout step.
- Interrupted transfer: run the same rsync command again; `--partial` allows a
  safer restart.
- Older rsync without `--info=progress2`: replace it with `--progress`.
- Never add `--delete` unless deletions have been previewed and explicitly
  approved; it can remove local-only results.
