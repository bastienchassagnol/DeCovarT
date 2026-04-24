#!/usr/bin/env bash

set -euo pipefail

NUMBER="${NUMBER:-}"
EVENT="${EVENT:-synchronize}"

if [[ -z "${NUMBER}" ]]; then
  echo "NUMBER environment variable is required."
  exit 1
fi

PR_JSON="$(gh pr view "${NUMBER}" --json title,body,author,labels,draft)"
AUTHOR="$(echo "${PR_JSON}" | jq -r '.author.login')"
TITLE="$(echo "${PR_JSON}" | jq -r '.title')"
BODY="$(echo "${PR_JSON}" | jq -r '.body // ""')"
LABEL_COUNT="$(echo "${PR_JSON}" | jq '.labels | length')"
DRAFT="$(echo "${PR_JSON}" | jq -r '.draft')"

# Skip bot-authored PRs.
if [[ "${AUTHOR}" == *"[bot]"* ]]; then
  echo "PR #${NUMBER}: bot author detected, skipping checks."
  exit 0
fi

REGEX_TITLE='^(feat|fix|chore|docs|refactor|revert|style|test|ci|build|perf|deps)(\(.*\))?(!)?: .+$'
REGEX_ISSUE='(?i)(close|closes|closed|fix|fixes|fixed|resolve|resolves|resolved)[[:space:]]+#[0-9]+'

ERRORS=""
WARNINGS=""

if [[ "${DRAFT}" != "true" && ! "${TITLE}" =~ ${REGEX_TITLE} ]]; then
  ERRORS+="- Pull request title should follow Conventional Commits, e.g. \`feat(scope): short description\`.\n"
fi

if [[ "${LABEL_COUNT}" -eq 0 ]]; then
  ERRORS+="- Please add at least one label to this pull request.\n"
fi

if [[ ("${EVENT}" == "opened" || "${EVENT}" == "ready_for_review") && ! "${BODY}" =~ ${REGEX_ISSUE} ]]; then
  WARNINGS+="- Recommended: reference an issue in the PR description (e.g. \`Closes #123\`).\n"
fi

if [[ -n "${ERRORS}${WARNINGS}" ]]; then
  {
    echo "@${AUTHOR} Thanks for this pull request."
    echo
    if [[ -n "${ERRORS}" ]]; then
      echo "**Required fixes:**"
      echo
      printf "%b" "${ERRORS}"
      echo
    fi
    if [[ -n "${WARNINGS}" ]]; then
      echo "**Recommended updates:**"
      echo
      printf "%b" "${WARNINGS}"
    fi
  } >pr_comment_body.txt

  gh pr comment "${NUMBER}" --body-file pr_comment_body.txt
  rm -f pr_comment_body.txt
fi

if [[ -n "${ERRORS}" ]]; then
  exit 1
fi

echo "PR #${NUMBER}: checks passed."
