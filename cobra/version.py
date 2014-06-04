#!/usr/bin/env python
from __future__ import print_function

__all__ = ("get_version")

from subprocess import check_output
from os import path, name, devnull
current_directory = path.dirname(path.abspath(__file__))
version_file = path.join(current_directory, "VERSION")

git_command = "git"
if name == "nt":
    def find_git_on_windows():
        from subprocess import CalledProcessError
        # first see if git is in the path
        try:
            check_output(["where", "/Q", "git"])
            # if this command succeeded, git is in the path
            return "git"
        # catch the exception thrown if git was not found
        except CalledProcessError:
            None
        # There are several locations git.exe may be hiding
        possible_locations = []
        from os import environ, listdir
        # look in program files for msysgit
        if "PROGRAMFILES(X86)" in environ:
            possible_locations.append("%s/Git/cmd/git.exe" % \
                environ["PROGRAMFILES(X86)"])
        if "PROGRAMFILES" in environ:
            possible_locations.append("%s/Git/cmd/git.exe" % \
                environ["PROGRAMFILES"])
        # look for the github version of git
        if "LOCALAPPDATA" in environ:
            github_dir = "%s/GitHub" % environ["LOCALAPPDATA"]
            if path.isdir(github_dir):
                for subdir in listdir(github_dir):
                    if not subdir.startswith("PortableGit"):
                        continue
                    possible_locations.append("%s/%s/bin/git.exe" % \
                        (github_dir, subdir))
        for possible_location in possible_locations:
            if path.isfile(possible_location):
                return possible_location
        # git was not found
        return "git"

    git_command = find_git_on_windows()


def call_git_describe(abbrev=7):
    try:
        with open(devnull, "w") as fnull:
            arguments = [git_command, "describe", "--tags",
                         "--abbrev=%d" % abbrev]
            return check_output(arguments, cwd=current_directory,
                                stderr=fnull).decode("ascii").strip()
    except:
        return None


def read_release_version():
    try:
        with open(version_file, "r") as infile:
            version = str(infile.read().strip())
        if len(version) == 0:
            version = None
        return version
    except:
        return None


def get_version():
    """Tracks the version number.

    The file VERSION holds the version information. If this is not a git
    repository, then it is reasonable to assume that the version is not
    being incremented and the version returned will be the release version as
    read from the file.

    However, if the script is located within an active git repository,
    git-describe is used to get the version information.

    The file VERSION will need to be changed by manually. This only
    needs to occur twice per release:
      - Once right before running git tag (set to the same as the version in
        the tag).
      - Once right after running git tag (set to next_version.dev)
    """

    git_version = call_git_describe()
    if git_version is None:  # not a git repository
        return read_release_version()
    return git_version


if __name__ == "__main__":
    print(get_version())
