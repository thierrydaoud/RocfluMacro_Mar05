# The following code is the most recent version of combined Hydrocode <u>RocfluMP</u> and particle-tracker <u>ppiclF</u>.

## 📘 **Documentations**: 
### RocfluMP
You can find the full RocfluMP manual [here](docs/Rocflu_manual.pdf).

### ppiclF
Documentation Website for ppiclF [here](https://dpzwick.github.io/ppiclF-doc/user/external.html).

## Getting Started:

## 🚀 Getting Started

Follow these steps to get a copy of the project up and running on your machine or HiPerGator account.

### 1. 🍴 Fork the Repository

- Click the **Fork** button at the top-right of this page to create your own copy of the repository.

### 2. 📥 Clone Your Fork

#### Option A: Clone via **HTTPS**
```bash
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name
```

#### Option B: Clone via **SSH** (recommended for contributors)

1. **Generate an SSH key** (if you don’t already have one):

   ```bash
   ssh-keygen -t ed25519 -C "your_email@example.com"
   ```

2. **Add your SSH key to the ssh-agent**:

   ```bash
   eval "$(ssh-agent -s)"
   ssh-add ~/.ssh/id_ed25519
   ```

3. **Add the SSH key to your GitHub account**:
   Display your public key with:
   ```bash
   cat ~/.ssh/id_ed25519.pub
   ```
   Go to GitHub → Settings → SSH and GPG keys → New SSH key
   
   Paste the copied key into the field and give it a descriptive title (e.g., My Laptop or HiPerGator).

4. **Clone your fork using SSH**:
   ```bash
   git clone git@github.com:your-username/your-repo-name.git
   cd your-repo-name
   ```

## 🛠️ Making Changes with Git

Once you’ve cloned your fork of the repository, here are the basic Git commands you’ll use when contributing to the project:

### 🔄 1. Pull the Latest Changes

Before making changes, make sure your local branch is up to date with the main project:

```bash
git fetch upstream
git checkout main
git merge upstream/main
```

### 📝 2. Create a New Branch for Your Changes
It's good practice to make changes on a separate branch, not directly on `main`.
```bash
git checkout -b feature/your-branch-name
```

### 💻 3. Make Your Code Changes
Edit the files you need using your favorite editor or IDE.

### ✅ 4. Stage and Commit Your Changes
```bash 
git add .
git commit -m "Descriptive message of what you changed"
```

### 📤 5. Push to Your Fork
```branch
git push origin feature/your-branch-name
```

### 🚀 6. Open a Pull Request
1. Go to your forked repository on GitHub.

2. You’ll see a prompt to open a pull request (PR) from your branch.

3. Add a clear description of your changes and submit.


## 🧼 Recommended Tips


- Use `.gitignore` (already found in your cloned directory) to avoid committing unwanted files.

- Run `git status` often to check what’s staged or modified.

- Run `git log --oneline` to view your commit history.
