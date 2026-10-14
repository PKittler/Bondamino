# Transparency: use of AI in this project

Bondamino was developed with the help of artificial intelligence. This file says how, who was
responsible for what, and how the results were checked.

## Summary

- Bondamino was created with **Claude**, an AI model by [Anthropic](https://www.anthropic.com/),
  used through **Claude Code**, Anthropic's AI coding assistant that reads, writes and tests code
  in a project folder.
- Claude wrote most of the source code, the automated tests, the translations and the
  documentation, including this file.
- A human, Philipp Kittler, directed the work: he defined the goals and features, made all
  design and licensing decisions, reviewed and tried out the results, and decided what was
  accepted into the project.
- Bondamino itself contains **no AI**. When you use Bondamino, your structures and data are
  processed only by Bondamino's own code on your server and are not sent to an AI service.

## How Bondamino was made

| Period | What happened |
|---|---|
| April 2025 | Project start with Django's standard project template (commits by the project's human authors). |
| From October 2026 | Development of the application with Claude Code: structure parsing, graph, reference values, statistics, Ramachandran evaluation, work groups, notes, exports, installation, security review, code comments and documentation. The model used was Claude Opus 5.5. |

Commits created with Claude are marked in the Git history with the line
`Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Division of responsibilities

**The maintainer (human)**

- set the scientific scope and the requirements for each feature,
- chose between the options Claude proposed (e.g. licence, databases, user interface),
- tried out the application and gave feedback on it,
- reviewed the changes and merged them into the project,
- is responsible for the project and its publication.

**Claude (AI)**

- wrote and changed the code (Python/Django, JavaScript, HTML/CSS, Docker),
- wrote the automated tests and ran them,
- researched reference data and licences and prepared them,
- translated the user interface and wrote the documentation,
- reviewed the code for security problems and outdated dependencies.

## How the results were checked

- **Automated tests:** about 120 tests check parsing, graph building, reference values, the
  Ramachandran evaluation, permissions, sign-in, invitations, configuration and setup. They
  are run before every change is accepted.
- **Comparison with established tools:** reference values come from published sources
  (Engh & Huber 2001, CCP4 Monomer Library, Top8000; see [CITATION.md](CITATION.md)). The
  Ramachandran evaluation follows the MolProbity method and was compared with the wwPDB
  validation results.
- **Real structures:** during development, every feature was checked with real PDB entries
  (e.g. 1LYZ, 1D3Z, 1EJG, 1ATP, 1HZH), including rendered pages and exported PDF files.
- **Security:** a security review (sign-in throttling, dependency audit with `pip-audit` and
  `npm audit`, safe output of file contents) was carried out; see the commit history.

## Limits

AI models can make mistakes that look plausible, both in code and in scientific details.
Despite the checks above, errors cannot be ruled out. Before you rely on a result in a
publication, check it for plausibility, e.g. by comparing with the wwPDB validation report of
the structure or with other validation software. If you find an error, please report it in
the project's issue tracker.

## Licence of AI-assisted content

The AI-assisted code and documentation are published under the same licence as the rest of
Bondamino; see [LICENSE.md](LICENSE.md).
