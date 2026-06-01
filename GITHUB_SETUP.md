# GitHub Setup Instructions

Your local Git repository is ready! Follow these steps to create a private repository on GitHub.

## ✓ Already Completed:

- ✓ Git repository initialized
- ✓ .gitignore configured (excludes BAM files, output folders, etc.)
- ✓ Initial commit created (13 files, 2533 lines)
- ✓ Git identity configured

## Next Steps: Push to GitHub

### Option 1: Using GitHub Website (Easiest)

1. **Go to GitHub:** https://github.com/new

2. **Create New Repository:**
   - Repository name: `methylation-decoder` (or your preferred name)
   - Description: `Automated pipeline for decoding binary messages from DNA methylation patterns`
   - **Privacy: ✓ Private** ← IMPORTANT
   - **DO NOT** initialize with README, .gitignore, or license (we already have these)
   - Click "Create repository"

3. **Push Your Code:**

   GitHub will show you commands. Use these:

   ```bash
   # In your methylation_decoder directory:
   
   git remote add origin https://github.com/YOUR_USERNAME/methylation-decoder.git
   git branch -M main
   git push -u origin main
   ```

   Replace `YOUR_USERNAME` with your GitHub username.

4. **Authentication:**
   
   GitHub will prompt for credentials. Use:
   - Username: Your GitHub username
   - Password: **Personal Access Token** (NOT your GitHub password)

   To create a token:
   - Go to: https://github.com/settings/tokens
   - Click "Generate new token" → "Generate new token (classic)"
   - Select scopes: ✓ `repo` (full control)
   - Click "Generate token"
   - Copy the token (you won't see it again!)
   - Use this as your password when pushing

### Option 2: Using GitHub CLI (gh)

If you have `gh` installed:

```bash
# Login to GitHub
gh auth login

# Create private repository and push
gh repo create methylation-decoder --private --source=. --push

# Done!
```

### Option 3: Using SSH (More secure, no tokens needed)

1. **Generate SSH key** (if you don't have one):

   ```bash
   ssh-keygen -t ed25519 -C "everythingblocked0101@gmail.com"
   # Press Enter for default location
   # Press Enter for no passphrase (or set one)
   ```

2. **Add SSH key to GitHub:**

   ```bash
   # Copy your public key
   cat ~/.ssh/id_ed25519.pub
   ```

   - Go to: https://github.com/settings/keys
   - Click "New SSH key"
   - Paste the key
   - Click "Add SSH key"

3. **Create repository on GitHub** (same as Option 1, step 2)

4. **Push using SSH:**

   ```bash
   git remote add origin git@github.com:YOUR_USERNAME/methylation-decoder.git
   git push -u origin main
   ```

## After Pushing:

Your private repository will be at:
```
https://github.com/YOUR_USERNAME/methylation-decoder
```

### Set Repository Details:

1. Go to your repository on GitHub
2. Click "Settings"
3. Add:
   - Description: "Automated pipeline for decoding DNA methylation patterns"
   - Topics: `bioinformatics`, `methylation`, `nanopore`, `dna-storage`, `python`
   - Website: (optional)

### Invite Collaborators (Optional):

For a private repo, to give access to others:
1. Go to Settings → Collaborators
2. Click "Add people"
3. Enter their GitHub username

## Repository Structure on GitHub:

```
methylation-decoder/
├── README.md                    # Main page
├── QUICKSTART.md               # Quick start guide  
├── LICENSE                     # MIT License
├── .gitignore                  # Configured
├── requirements.txt            # Python dependencies
├── master_pipeline.py          # Main entry point
├── example_run.sh              # Example script
├── lib/                        # Libraries
│   ├── dynamic_threshold_detector.py
│   └── brick_position_extractor.py
└── scripts/                    # Pipeline steps
    ├── index_bam.py
    ├── generate_bed.py
    ├── detect_threshold.py
    └── decode_sequence.py
```

## Future Updates:

When you make changes:

```bash
# Check what changed
git status

# Add changes
git add .

# Commit with message
git commit -m "Description of changes"

# Push to GitHub
git push
```

## Troubleshooting:

### "Permission denied (publickey)"
- You need to set up SSH key (see Option 3 above)
- Or use HTTPS with personal access token (Option 1)

### "remote: Repository not found"
- Check repository name matches
- Ensure you're using correct GitHub username
- For private repos, make sure you have access

### "Authentication failed"
- For HTTPS: Use **personal access token**, not password
- For SSH: Make sure key is added to GitHub

### "refused to merge unrelated histories"
```bash
git pull origin main --allow-unrelated-histories
# Then push again
```

## Clone to Other Machines:

After pushing, clone on other computers:

```bash
# HTTPS (requires token)
git clone https://github.com/YOUR_USERNAME/methylation-decoder.git

# SSH (requires SSH key)
git clone git@github.com:YOUR_USERNAME/methylation-decoder.git
```

## Keep Repository Private:

To verify it's private:
1. Go to repository on GitHub
2. Check for 🔒 **Private** badge next to repository name
3. Settings → Danger Zone → Change visibility (if needed)

## Recommended: Add Repository Badges

Add to top of README.md:

```markdown
![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Status](https://img.shields.io/badge/status-production-brightgreen.svg)
```

---

**Your repository is ready to push!** 🚀

Choose Option 1 (Website + HTTPS) if you're new to GitHub.
Choose Option 3 (SSH) for long-term convenience.
