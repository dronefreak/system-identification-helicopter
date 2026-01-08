# Security Policy

## Supported Versions

We currently support the following versions of this project with security updates:

| Version | Supported          |
| ------- | ------------------ |
| Latest  | :white_check_mark: |
| < 1.0   | :x:                |

## Reporting a Vulnerability

We take the security of our system identification toolkit seriously. If you discover a security vulnerability, please follow these guidelines:

### How to Report

**Please do NOT report security vulnerabilities through public GitHub issues.**

Instead, please report them via one of the following methods:

1. **Email**: Send details to the project maintainers (check CONTRIBUTORS.md for contact information)
2. **GitHub Security Advisories**: Use the [Security Advisories](../../security/advisories/new) feature on GitHub (recommended)

### What to Include

When reporting a vulnerability, please include:

- **Description**: Clear description of the vulnerability
- **Impact**: Potential impact and attack scenarios
- **Reproduction**: Detailed steps to reproduce the issue
- **MATLAB Version**: Version of MATLAB where the issue was found
- **Operating System**: OS and version information
- **Affected Files**: Specific files or functions affected
- **Suggested Fix**: If you have a proposed solution (optional)

### Example Report

```
Subject: [SECURITY] Potential Command Injection in data loading function

Description:
The load_flight_data.m function may be vulnerable to command injection
when processing user-supplied file paths.

Impact:
An attacker could potentially execute arbitrary MATLAB code by crafting
a malicious file path.

Steps to Reproduce:
1. Call load_flight_data() with a specially crafted path
2. [Additional steps]

MATLAB Version: R2020a
Operating System: Windows 10

Affected Files:
- src/utils/load_flight_data.m (line 45)

Suggested Fix:
Implement strict input validation using validateattributes() and
sanitize file paths before passing to load().
```

### Response Timeline

- **Initial Response**: Within 48 hours
- **Assessment**: Within 1 week
- **Fix Timeline**: Depends on severity
  - **Critical**: Within 7 days
  - **High**: Within 30 days
  - **Medium**: Within 90 days
  - **Low**: Next scheduled release

### Security Update Process

Once a vulnerability is confirmed:

1. We will acknowledge receipt and begin investigation
2. We will work on a fix and keep you informed of progress
3. We will prepare a security advisory
4. We will release a patched version
5. We will publicly disclose the vulnerability (after fix is available)

### Disclosure Policy

- **Coordinated Disclosure**: We follow responsible disclosure practices
- **Credit**: We will credit researchers who report vulnerabilities (unless anonymity is requested)
- **Public Disclosure**: After a fix is released and users have had time to update

## Security Best Practices

When using this toolkit:

### Data Security

- **Never commit sensitive data**: Use `.gitignore` for flight data containing proprietary information
- **Use Git LFS**: For large data files (see docs/GIT_LFS_SETUP.md)
- **Sanitize outputs**: Review generated reports for sensitive information before sharing

### Code Execution

- **Validate inputs**: Always validate user inputs and file paths
- **Sandbox execution**: Run optimization algorithms in isolated environments when processing untrusted data
- **Review configurations**: Check configuration files before running experiments

### File Operations

- **Path validation**: Use absolute paths and validate file locations
- **Permission checks**: Ensure proper file permissions for data directories
- **Temporary files**: Clean up temporary files after execution

### MATLAB-Specific

- **Eval safety**: This toolkit does not use `eval()` or `evalin()` with user input
- **Load safety**: Data loading functions validate file types and structure
- **Parallel computing**: When using parallel pools, ensure worker environments are secure

## Known Security Considerations

### Current Limitations

1. **File Path Handling**: File paths are validated but may need additional sanitization on some platforms
2. **Data Validation**: Flight data validation checks structure but not content authenticity
3. **Configuration Files**: Configuration files are trusted - only load configs from trusted sources

### Recommendations

- Run MATLAB with appropriate user permissions (avoid running as administrator/root)
- Keep MATLAB and required toolboxes updated
- Review third-party dependencies for security updates
- Use version control best practices (see .gitattributes and .gitignore)

## Security Enhancements Roadmap

Future security improvements planned:

- [ ] Digital signatures for configuration files
- [ ] Checksum validation for data files
- [ ] Encrypted storage options for sensitive parameters
- [ ] Audit logging for experiment execution
- [ ] Integration with MATLAB security features

## Compliance and Standards

This project follows:

- **MATLAB Coding Standards**: Secure coding practices for MATLAB
- **OWASP Guidelines**: Where applicable to scientific computing
- **Academic Integrity**: Reproducibility and data integrity standards

## Third-Party Dependencies

This project uses:

- **MATLAB**: R2014b or later (keep updated for security patches)
- **Toolboxes**: Parallel Computing Toolbox, Optimization Toolbox (optional)
- **Git LFS**: For large file management (keep updated)

Monitor these dependencies for security updates.

## Contact

For security-related questions or concerns:

- Review this document first
- Check existing security advisories
- Contact maintainers privately for new issues
- Use GitHub Discussions for general security questions (not vulnerabilities)

---

**Last Updated**: 2026-01-08

Thank you for helping keep this project and its users safe!
