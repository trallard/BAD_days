---
layout: page
permalink: /3-outline/
title: Day 3
category: module
---
**Prerequisites:**

Addressed to researchers already familiar with R and basics of statistics

The following packages are required:

- RNAseq123
- edgedR
- Mus.musculus
- limma
- Glimma

....plus plotting packages

**Description:**

In the BAD DAY 3 we will present the Bioconductor work-flow R package RNAseq123 from <http://www.bioconductor.org/help/workflows/RNAseq123/>.

We will in particular discuss how to:

- Import and organise data
- Preprocess data
- Perform exploratory plots
- Analyse Differential Expression in genes
- Perform Exploratory analysis and graphical representation of DE genes
- Test for Gene set enrichment

---

## Outline/ content


<ul >
{% assign pages_list = site.pages %}
{% for node in pages_list %}
{% if node.title != null %}
{% if node.layout == "default" %}
{%if node.tags %}

{% for tag in node.tags %}
{% if tag == 'Day3' %}
<!-- Note you need to prepend the site.baseurl always-->
<li><a href="{{ node.url | absolute_url}}">{{ node.title }}</a>
</li>
{% endif %}
{% endfor %}

{% endif %}
{%endif%}
{% endif %}
{% endfor %}
</ul>




<a href="{{site.url}}{{site.baseurl}}/index.html" class="float" >
<i class="fa fa-home my-float"></i>
</a>
