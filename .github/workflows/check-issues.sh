#!/usr/bin/env bash

set -euo pipefail

NUMBER="${NUMBER:-}"
if [[ -z "${NUMBER}" ]]; then
  echo "NUMBER environment variable is required."
  exit 1
fi

ISSUE_JSON="$(gh issue view "${NUMBER}" --json number,title,author,labels)"
AUTHOR="$(echo "${ISSUE_JSON}" | jq -r '.author.login')"
LABEL_COUNT="$(echo "${ISSUE_JSON}" | jq '.labels | length')"

WARNINGS=""

if [[ "${LABEL_COUNT}" -eq 0 ]]; then
  WARNINGS+="- Please add at least one label to this issue.\n"
fi

if [[ -n "${WARNINGS}" ]]; then
  {
    echo "@${AUTHOR} Thanks for opening this issue."
    echo
    echo "Please update it with the following:"
    echo
    printf "%b" "${WARNINGS}"
  } >issue_comment_body.txt

  gh issue comment "${NUMBER}" --body-file issue_comment_body.txt
  rm -f issue_comment_body.txt
  exit 1
fi

echo "Issue #${NUMBER}: checks passed."
