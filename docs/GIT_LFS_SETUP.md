# Git LFS Setup Guide

This guide explains how to set up Git Large File Storage (LFS) for managing large `.mat` data files in the repository.

## Why Git LFS?

Git LFS is designed for versioning large files efficiently. Instead of storing large binary files directly in the Git repository, LFS stores pointers to the actual files, which are stored separately.

**Benefits**:
- Faster cloning and fetching
- Reduced repository size
- Better performance for large binary files
- Efficient storage of file versions

## Current Data Files

The repository contains `.mat` files ranging from 1.5KB to 30MB:
- Small files (< 1MB): Can be tracked normally
- Medium files (1-5MB): Recommended for LFS
- Large files (> 5MB): **Should use LFS**

## Installation

### 1. Install Git LFS

**macOS** (using Homebrew):
```bash
brew install git-lfs
```

**Windows**:
- Download from https://git-lfs.github.com/
- Run the installer

**Linux** (Ubuntu/Debian):
```bash
sudo apt-get install git-lfs
```

**Linux** (Fedora/CentOS):
```bash
sudo yum install git-lfs
```

### 2. Initialize Git LFS

In your repository:
```bash
git lfs install
```

This only needs to be done once per machine.

## Configuring LFS for This Repository

The repository is already configured to track `.mat` files with Git LFS via `.gitattributes`:

```
*.mat filter=lfs diff=lfs merge=lfs -text
```

However, existing files are **not automatically** converted to LFS. See migration section below.

## For New Contributors

### Clone with LFS

```bash
# Clone repository
git clone https://github.com/dronefreak/system-identification-helicopter.git
cd system-identification-helicopter

# LFS files are automatically downloaded
# If not, manually fetch them:
git lfs pull
```

### Check LFS Status

```bash
# See which files are tracked by LFS
git lfs ls-files

# Check LFS status
git lfs status
```

## Migrating Existing Files to LFS

If you want to convert existing `.mat` files to LFS tracking:

⚠️ **Warning**: This rewrites Git history. Coordinate with all contributors first.

```bash
# Navigate to repository root
cd /path/to/system-identification-helicopter

# Migrate existing .mat files to LFS
git lfs migrate import --include="*.mat" --everything

# Force push (required after history rewrite)
git push --force origin main
```

### Selective Migration

To migrate only large files (e.g., > 5MB):

```bash
# Find large .mat files
find data/ -name "*.mat" -size +5M

# Migrate specific files
git lfs migrate import \
  --include="data/experiments/11_05_09.bin-1077528.mat,data/experiments/SysID405705308.mat" \
  --everything
```

## Working with LFS Files

### Adding New Files

New `.mat` files are automatically tracked by LFS (if `.gitattributes` is configured):

```bash
# Add new .mat file
git add data/experiments/new_data.mat
git commit -m "Add new experiment data"
git push
```

### Checking File Status

```bash
# List files tracked by LFS
git lfs ls-files

# Example output:
# 4c7e2a8f * data/experiments/best.mat
# 9f3b1d2e * data/experiments/best2.mat
```

### Fetching LFS Files

```bash
# Fetch all LFS files
git lfs fetch

# Fetch and checkout LFS files
git lfs pull

# Fetch LFS files for specific remote
git lfs fetch origin
```

## Storage and Bandwidth

### LFS Storage Limits

- GitHub Free: 1 GB storage, 1 GB/month bandwidth
- GitHub Pro: 2 GB storage, 2 GB/month bandwidth
- Additional packs available for purchase

### Checking Usage

```bash
# On GitHub, check: Settings > Billing > Git LFS Data
```

### Optimizing Storage

1. **Remove old LFS files** from history if no longer needed
2. **Use data compression** before committing large files
3. **Archive old experiments** to external storage
4. **Prune old LFS objects**:
   ```bash
   git lfs prune
   ```

## Alternative: External Data Storage

For very large datasets or limited LFS storage, consider:

1. **Cloud Storage**: S3, Google Cloud Storage, Azure Blob
2. **Institutional Storage**: University/research institution servers
3. **Data Repositories**: Zenodo, Figshare, Dryad
4. **Version Control**: DVC (Data Version Control)

Store data externally and keep only:
- Small sample/test datasets in repository
- Download scripts for full datasets
- Documentation of data sources

## Troubleshooting

### LFS Files Not Downloaded

```bash
# Check LFS configuration
git lfs env

# Pull LFS files
git lfs pull
```

### Large Clone Size

LFS files may not be the issue. Check:
```bash
# Find large non-LFS files
git rev-list --objects --all | \
  git cat-file --batch-check='%(objecttype) %(objectname) %(objectsize) %(rest)' | \
  sed -n 's/^blob //p' | \
  sort --numeric-sort --key=2 | \
  tail -20
```

### Converting Files Back from LFS

```bash
# Untrack .mat files from LFS
git lfs untrack "*.mat"

# Migrate back to regular Git
git lfs migrate export --include="*.mat"
```

## Best Practices

1. **Track large files with LFS** (> 1MB recommended)
2. **Don't track frequently changing large files** (LFS stores all versions)
3. **Use compression** for large data when appropriate
4. **Document data sources** and preprocessing steps
5. **Use semantic versioning** for data releases
6. **Add checksums** for important datasets
7. **Test LFS setup** before pushing large files
8. **Communicate with team** before migrating to LFS

## Current Status

**Repository Status**: `.mat` files configured for LFS tracking in `.gitattributes`

**Migration Status**: Existing files NOT yet migrated (requires history rewrite)

**Recommendation**:
- For new development: LFS will work automatically
- For existing files: Consider migrating files > 5MB
- Coordinate with all contributors before migrating

## See Also

- Official Git LFS Documentation: https://git-lfs.github.com/
- GitHub LFS Guide: https://docs.github.com/en/repositories/working-with-files/managing-large-files
- `docs/DATA_CATALOG.md` - Data file documentation
- `src/utils/load_data.m` - Data loading utilities
