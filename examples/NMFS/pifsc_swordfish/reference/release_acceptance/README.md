# Local swordfish release acceptance

The clean isolated run `runs/swordfish-release-003` passed all ten stages.
`local_acceptance_summary.json` records actual compiler/platform, stage times,
manifest identity and the remaining scientific review boundary. All 116 build
artifact hashes and source identities were rechecked after publication.

The workspace runtime/model source and the frozen accepted source match by
content hash. All eight original inputs retain their recorded hashes. Local
portable AddressSanitizer/UndefinedBehaviorSanitizer checks passed; their log
is retained here. Hosted Ubuntu/macOS CI is configured but was not run in this
session. Failed development snapshots remain unpublished under ignored runs/.

Documentation: `documentation_validation.json` binds the final 41-page PDF and
280-declaration API catalog to their hashes. All pages were rendered/reviewed.
See `../../docs/RELEASE_READINESS.md` for operational guarantees and limitations.
