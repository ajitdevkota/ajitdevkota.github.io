---
title: Developing Asset Maintenance Application
date: 2021-01-31 16:41:58 +/-0000
categories: [Blogging, Demonstration]
tags: [xamarin, c#, android, ios, mysql]
---

## Overview
During this project, I worked in a team to develop an application for a company called ENGIE,
which looked on ways to simplify their assets maintaining process. Furthermore, the company
wanted to reduce its Carbon Footprint by going entirely paperless. The client required us to develop an application
that supports both Android & iOS devices.

![Window shadow](/assets/img/posts/engie-home.png){: .shadow style="max-width: 50%" }
_App's home page_

**Workflows used**: Git and Waterfall.
**Technologies used**: Xamarin, C#, XAML, AWS, MySQL.

**Note**: *This project was done as part of assessment for module CSC2033.*

## Customer Requirements
Build an App that ENGIE’s subcontractors can use that will:

1. Link with QR code or barcode to identify the specific asset that is subject to a service (QR/barcode will be linked to asset register).
2. Depending on the asset, serve up the corresponding compliance form.
3. In all cases serve up a asset condition survey form.
4. Compliance with Data Protection Act 2018.
5. Comply with ENGIE branding guidelines.

## Product Development
We decided to use the Xamarin mobile app platform to create a cross-platform mobile app that shares a native user
interface and most code between different mobile operating systems. All coding was written using C# programming language.
The Graphical User Interface (GUI) of the application was made using XAML user interface mark-up language.

## My Roles and Responsibilities
- Database connection
- Creation of Login & Registration pages
- Biometric authentication implementation (Fingerprint & Face ID)
- Security/Encryption
- PDF support implementation
- GUI Creation & Design (XAML)
- Dark Mode support
- Testing

## Learning Outcomes
- Designed and implemented a large software system.
- Gained practical experience in real-world software development and team working.
- Demonstrated practical experience in issues such as team structure, document preparation, project management.
- Fulfilled appropriate roles within a team and applied particular skills to the job in hand.
- Gained an experience of working with appropriate industrial-strength tool-chains.
- Developed employability skills.

## Learn More
For more knowledge about the project, visit my project's [repository](https://github.com/mantas2000/Engie-Asset-Maintenance-App){: target="_blank"}.
