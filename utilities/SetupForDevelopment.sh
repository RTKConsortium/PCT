#!/usr/bin/env bash
#==========================================================================
#
#   Copyright NumFOCUS
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#          https://www.apache.org/licenses/LICENSE-2.0.txt
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
#==========================================================================*/


# Run this script to set up the git hooks for committing changes to PCT.
# For more information, see:
#   https://www.itk.org/Wiki/ITK/Git#Hooks
#   https://www.itk.org/Wiki/Git/Hooks

#!/usr/bin/env bash
#=============================================================================
# Copyright 2010-2012 Kitware, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#=============================================================================

# Run this script to set up local Git hooks for this project.

# Project configuration instructions:
#
# - Publish a "hooks" branch in the project repository such that
#   clones will have "refs/remotes/origin/hooks".
#
# - Populate adjacent "config" file with:
#    hooks.url    = Repository URL publishing "hooks" branch
#    hooks.branch = Repository branch instead of "hooks"

egrep-q() {
	egrep "$@" >/dev/null 2>/dev/null
}

die() {
	echo 1>&2 "$@" ; exit 1
}

# Make sure we are inside the repository.
cd "${BASH_SOURCE%/*}" &&

# Select a hooks branch.
if url=$(git config --get hooks.url); then
	# Fetch hooks from locally configured repository.
	branch=$(git config hooks.branch || echo hooks)
elif git for-each-ref refs/remotes/origin/hooks 2>/dev/null |
     egrep-q 'refs/remotes/origin/hooks$'; then
	# Use hooks cloned from origin.
	url=.. && branch=remotes/origin/hooks
elif url=$(git config -f config --get hooks.url); then
	# Fetch hooks from project-configured repository.
	branch=$(git config -f config hooks.branch || echo hooks)
elif git ls-remote 'https://github.com/RTKConsortium/PCT' 2>/dev/null |
	egrep -q 'refs/heads/hooks$'; then
	# Fallback on PCT's hooks branch
	url='https://github.com/RTKConsortium/PCT' && branch='refs/heads/hooks'
else
	die 'This project is not configured to install local hooks.'
fi &&

echo $url
echo $branch

# Populate ".git/hooks".
echo 'Setting up git hooks...' &&
git_dir=$(git rev-parse --git-dir) &&
mkdir -p "$git_dir/hooks" &&
cd "$git_dir/hooks" &&
if ! test -e .git; then
	git init -q || die 'Could not run git init for hooks.'
fi &&
git fetch -q "$url" "$branch" &&
git reset -q --hard FETCH_HEAD || die 'Failed to install hooks'

cd ../..

# Set up KWStyle hook.
echo "Setting up the KWStyle hook..."
git config hooks.KWStyle.conf "cmake/KWStyle/PCT.kws.xml"
git config hooks.KWStyle.overwriteRulesConf "cmake/KWStyle/PCTOverwrite.txt"

echo "Done."
