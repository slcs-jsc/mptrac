# Release management


MPTRAC and JURASSIC follow a coordinated, lightweight release process.
Both projects aim for approximately two regular releases per year,
preferably in January and July, when development activity is relatively low.
Release dates are flexible; the projects may also publish additional
releases after major development or validation milestones, such as a
comprehensive code review. The two projects should normally be released
around the same time, but a release of one project does not depend on the
readiness of the other.

Version numbers reflect the changes made, not the release date. Existing
Git tags use the historical `vX.Y` convention. These two-component tags
are retained for compatibility with previous releases; they are not full
three-component Semantic Versioning (`MAJOR.MINOR.PATCH`) identifiers.
For new versions, choose the appropriate change level based on compatibility
and functionality, and document the versioning convention consistently.
Do not rename historical tags.

A release normally needs only a short maintenance window of approximately
one or two days. Keep this focused on small fixes and cleanup rather than
introducing substantial new features.

### Release checklist

1. Review outstanding issues, resolve suitable small bugs, and record any
   known limitations that remain open.
2. Run the project's regression tests and check the GitHub Actions results.
   Run additional scientific or HPC validation when the changes warrant it.
3. Check documentation, citation metadata, and release-relevant version
   references. Confirm that the intended commit is ready for publication.
4. Select a version and create an annotated Git tag on the reviewed commit,
   for example `git tag -a vX.Y -m "Release vX.Y"`. Push the specific tag
   with `git push origin vX.Y`.
5. Summarize changes since the previous tag with
   `git log <previous-tag>..vX.Y --oneline`. Publish a GitHub Release
   from that tag with concise release notes and known limitations.
6. Check that Zenodo has archived the new release and provides a
   version-specific DOI. Use the project's concept DOI for references to
   the software across versions.
7. Register the new release in the JuSER publication database if applicable.

Release publication is a manual decision. A scheduled release window does
not require publishing code that is not ready.

