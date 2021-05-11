## Introduction to bioinformatics : Practice day 1

<div class="tables-start"></div>
|**Title**| Basic linux commands|
|---------|-------------------------------------------|
|**Training dataset:**| None |
|**Questions:**| <ul><li>How do I use the command line?</li><li>How do I navigate the file system?</li><li>How do I know which software I am using?</li><li>How do I use Nextflow with a Singularity image?</li></ul>|
|**Objectives**:|<ul><li>Remember how to use the command line.</li><li>Learn how to execute the Nextflow pipeline over our Singularity image for the following exercises.</li></ul>|
|**Time estimation**:| 30 min |
|**Key points**:|<ul><li>Remeber the shell basics</li><li>Learn how to call the pipeline for future exercises.</li></ul>|
<div class="tables-end"></div>

#### How do I use the command line?

Open a terminal and type into it. Please, remember to:
* use TAB to autocomplete and suggest commands, paths and files!
* use the left-right direction arrows to move the cursor inside the terminal!
* use the up-down direction arrows to navigate the command history!
* DON'T use white-spaces or specials characters (as Ñ) to name your files!

#### How do I know who and where am I?

It is important to know basic things as who is using the machine, which user you are connected as, which permissions you have and which OS is running:
```
who
whoami
id
uname -a
```

#### When am I?

Show today's date in the calendar, this year's december month and today's date and time:
```
cal
cal -m 7 2019
date
```

#### Clean the screen

Clean the terminal to work in an empty screen again:
```
clear
```

#### How do I navigate the file system?

Let's remember the basics: *pwd cd ls mkdir mv rm rmdir less nano*. We are going to use those commands to:

Check our working directory:

```
pwd
```

Move to our Desktop folder (using the absolute path and the alias "~", which means "path to your home folder"):

```
cd ~/Desktop
```

Create a directory called "myDir" (Linux is case sensitive and does not like white spaces in names):

```
mkdir myDir
```

Move to the new folder (using a relative pathway):

```
cd myDir
```

Check our working directory (always do it before executing something):

```
pwd
```

Create the folders "asdf", "AsDf", "ASDF" and "tmp" (at once, commands change their behavior depending on the parameters):

```
mkdir asdf AsDf ASDF tmp
```

Create a file inside "tmp" called "myFile.txt" (using a relative pathway, you can work with files outside your working directory):

```
nano tmp/myFile.txt
```

and write whatever you want and save it with __Ctrl + O__ and close the new file with __Ctrl + X__

Rename "myFile.txt" to "whateverIwant" (Linux does not require file extensions):

```
mv tmp/myFile.txt tmp/whateverIwant
```

See the contents of "whateverIwant":

```
less tmp/whateverIwant
```

Copy "whateverIwant" to a new file "whateverYOUwant", edit its contents with nano:

```
cd tmp
cp whateverIwant whateverYOUwant
nano whateverYOUwant
less whateverYOUwant
```

Concatenate both files into one:

```
cat whateverIwant whateverYOUwant
cat whateverIwant whateverYOUwant > whateverWEwant
less whateverWEwant
```

List all the files with their properties:

```
ls -lah
```

Remove the file "whateverIwant" (".." means parent directory):

```
cd ..
rm tmp/whateverIwant
```

Remove the folders inside "myDir" (wildcard character "\*" means "any character once or more times, or nothing"):

```
pwd
cd ..
rmdir ./*
rm myDir/tmp/whateverYOUwant myDir/tmp/whateverYOUwant
rmdir ./*
```

Go back to Desktop and remove everything you created (".." means parente directory, while "." refers to the directory itself):

```
cd ..
rm -rf myDir
```

Return to your home directory (without parameters, the behavior of the command changes):

```
cd
```

#### How do I know which software I am using?

Software may (and will) be installed in many different places. To discover the one you have loaded in your PATH use `which`, to see all the places where the shell is looking for software check the variable `$PATH`, to know the version of the software use the apropiate parameter (`-h --help -v --version`) and to check the manual of the software use `man`.

```
which git
echo $PATH
echo $USER
git -h
git --version
man git
```

#### How do I use Nextflow with a Singularity image?

```
nextflow
```

So, what now? In order to execute a nextflow pipeline, we need to tell it to `run` a project which contains a `main.nf` script written in groovy + the pipeline languages:

```
nextflow run /home/$USER/Documents/wgs/introduction_to_bioinformatics
```

Optionally, we can pass a config file, and specify the .nf script inside a project:

```
nextflow -C /home/$USER/Documents/wgs/bacterial_wgs_training/nextflow.config \
run /home/$USER/Documents/wgs/introduction_to_bioinformatics/main.nf
```

There is no need to download the software you want to execute, you can also execute a github repository:

```
nextflow run BU-ISCIII/introduction_to_bioinformatics
```

This is how we will execute the exercises during this course, so let's remove the downloaded repository to fre some space:

```
rm -rf /home/$USER/Documents/wgs/introduction_to_bioinformatics
```

Finally, let's ask how to use the pipeline:

```
nextflow run BU-ISCIII/introduction_to_bioinformatics --help
```

There is one big detail left. The software needed to execute the pipeline is no installed in our machine. Thankfully, we have a singularity image (container) ready for this course, and our pipeline has already being configurated to know where to find it and how to use it. Use the right argument and go for it:

```
nextflow run BU-ISCIII/introduction_to_bioinformatics -profile singularity
```
