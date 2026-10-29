"""Create a versioned draft with inherited files; never upload, delete or publish."""
import json
import os
from pathlib import Path
import re
import urllib.request
from datetime import datetime, timezone


def main():
    token = os.environ['ZENODO_TOKEN']
    source = os.environ['SOURCE_RECORD']
    version = os.environ['NEW_VERSION']
    if not source.isdigit() or not re.fullmatch(r'\d+\.\d+\.\d+', version):
        raise ValueError('Expected numeric source record ID and semantic version')

    def request(url, method='GET', payload=None):
        if not url.startswith('https://zenodo.org/api/deposit/depositions/'):
            raise ValueError('Unexpected deposition endpoint')
        headers = {'Authorization': 'Bearer ' + token, 'User-Agent': 'LayeredRock versioned archive',
                   'Accept': 'application/json', 'Content-Type': 'application/json'}
        data = None if payload is None else json.dumps(payload).encode()
        req = urllib.request.Request(url, data=data, method=method, headers=headers)
        with urllib.request.urlopen(req, timeout=120) as response:
            return json.load(response)

    original_url = 'https://zenodo.org/api/deposit/depositions/' + source
    original = request(original_url)
    if not original.get('submitted'):
        raise RuntimeError('Source must be an already published version')
    result = request(original_url + '/actions/newversion', method='POST')
    draft_url = result['links']['latest_draft']
    draft = request(draft_url)
    if draft.get('submitted') or str(draft['id']) == source:
        raise RuntimeError('Expected a distinct unpublished version draft')
    metadata = dict(draft['metadata'])
    metadata['version'] = version
    description = re.sub(r'GPL-3\.0(?!-only)', 'GPL-3.0-only', metadata.get('description', ''))
    revision_note = ('<p>Version 1.1.0 synchronizes the manuscript figure numbering while retaining '
        'the original data and script paths. It includes the 80 original prospective prediction JSON '
        'records (48 layout states and 32 new-length states), their saved method and batch records, '
        'and selected numerical clarification tables. The thirteen previously archived representative '
        'MPH models are retained. Source code is GPL-3.0-only; numerical data, prediction records and '
        'COMSOL MPH model files are CC BY 4.0. Original CT images retain their source attribution and '
        'license conditions. COMSOL and MATLAB software are not distributed in this archive.</p>')
    if 'Version 1.1.0 synchronizes the manuscript figure numbering' not in description:
        description += revision_note
    metadata['description'] = description
    draft = request(draft_url, method='PUT', payload={'metadata': metadata})
    doi = draft.get('metadata', {}).get('prereserve_doi', {}).get('doi') or draft.get('doi')
    receipt = {
        'prepared_utc': datetime.now(timezone.utc).isoformat(),
        'source_record_id': int(source), 'source_version': original['metadata'].get('version'),
        'new_deposition_id': draft['id'], 'new_version': version, 'reserved_doi': doi,
        'draft_api_url': draft_url,
        'draft_html_url': draft.get('links', {}).get('html'),
        'submitted': draft.get('submitted', False), 'final_publication_performed': False,
        'files_inherited': [{'name': f['filename'], 'bytes': f['filesize']} for f in draft.get('files', [])],
        'old_files_changed': False, 'new_files_uploaded': False,
        'record_level_license': draft['metadata'].get('license'),
        'description': draft['metadata'].get('description'),
    }
    Path('zenodo_newversion_draft_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2), flush=True)
    if not doi:
        raise RuntimeError('Draft created but no reserved DOI returned; inspect draft before packaging')


if __name__ == '__main__':
    main()
