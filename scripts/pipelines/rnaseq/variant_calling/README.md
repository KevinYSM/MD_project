/# Variant Calling Workflow
<!-- 
NF-core's RNASeq pipeline takes raw fastq files and produces filtered mapped reads files (*.bam, *.bai). The next step from this is variant calling.
There is a NF-core rnavar pipeline designed for this purpose, but I'll be doing it myself as:
1. The NF-core rnavar pipeline takes fastq files as input, requiring recomputation of steps already performed by NF-core RNASeq
2. I want to perform variant black listing as per a PDX paper I read (will cite here at a later date)
3. Good practice for me to write some *quality* pipelines 
4. The rnavar pipeline appears to be quite out of date, with its last release being in Jun 20, 2022. 
-->



<!--
Summary of Tools and Versions Used in the Pipeline:
GATK
Tabix (htslib?)
SnpEff
Ensembl VEP
 -->

<!--
link to your reader to your repository's bug page, and let them know if you're open to contributions
-->

## How to use Project Name:
<!--
Link to a webpage, web shell (e.g. runkit), or downloadable executable that demonstrates the project. 
        * note that when the reader is modifying the code, they will compare their modified version to the demo to see if their changes worked as they expected them to. Your demo is their reference
-->

### Installation:
<!--Explain how to import the modules of the project into the reader's codebase, install the containers of the project in the reader's cluster, or flash the binary of the project onto the reader's hardware-->

### API Methods | Modules:
<!--
List the methods or modules your project provides.
-->

## How Project Name works:
<!--
Explain how execution works. What is the entry point for your code? Which files correspond to which functionality? What is the lifecycle of your project? Are there any singletons, side effects or shared state among instances of your project? Take extra care to explain design decisions. After all, you wrote an ENTIRE codebase around your opinions. Make sure that the people using it understand them.
-->

## Roadmap:
<!--
List the releases that you have added to each project, and any future releases you would like to do. If there is a date for future release, put it here. If not, let people know that there is no defined timeframe for future releases.
-->

## Contribute to Project Name:
<!--
What are the prerequisites for contributing to the code?
    * provide users with containerized development environments, virtual machines, or, if developing for an embedded system, a pre-built OS image. Don't make them set up an environment from scratch.
-->

### Develop:
<!--
Tell your reader how to run the code in the development environment
-->

#### Repository Structure:
<!--
List each file, and what it does.
    * Identify whether you are open to pull requests for a specific file or not.
-->

| File or Folder | What does it do? | When should you modify it? |
| :------------- | :--------------- | :------------------------- |
|                |                  |                            |

### Test:
<!--
When the reader runs the code, what are the expected inputs and outputs?
How can the reader tell if the code is malfunctioning?
-->

### Document:
<!--
How should the reader document changes and additions to the code?
-->

### Deploy:
<!--
How is the code deployed? When the reader submits a pull request, how is the code merged into main and converted into a package?
-->

<!--
Additional tip: SHOW, don't TELL
* DON'T try to sell your reader on using your code. Don't spend words on clever analogies or context. That material is great for a blog post or video, but bad for the documentation included in repository. Your reader wants to run the code, not read about it. Help your reader get to 'hello world' as fast as possible.
* DO make diagrams. A diagram can help yoru reader organize information in ways that words alone can't.
    * Do not put more than 50 nodes and edges into a single diagram. It will turn into an indecipherable spaghetti-string mess. Keep diagrams simple.
-->