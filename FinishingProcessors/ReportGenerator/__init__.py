from jinja2 import Template
from Processor import Processor, ProcessorReturn


class ReportGenerator(Processor):
    """Generate a report for the dataset"""
    name = "ReportGenerator"
    help = "Generate a html report for the dataset"

    datasets = []

    def _process(self):
        template = Template("""
<html>
<head>
    <title>{{ bssFile.code }}</title>
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css" integrity="sha384-Vkoo8x4CGsO3+Hhxv8T/Q5PaXtkKtu6ug5TOeNV6gBiFeWPGFN9MuhOf23Q9Ifjh" crossorigin="anonymous">
    <style>
    .content {
  padding: 3rem 1.5rem;

}
    </style>
</head>
<body>
 <nav class="navbar navbar-expand-md navbar-dark bg-dark fixed-top">
  <a class="navbar-brand" href="#">{{ bssFile.code }}</a>
</nav>

<main role="main" class="container">
<div class="content">
{% if bssFile.depiction %}
<img src="{{bssFile.depiction}}" alt="Depiction of the object"/>
{% endif %}
<h2>Names</h2>
<ul>    
    {% for name in bssFile.names %}
        <li>{{ name }}</li>
    {% endfor %}
</ul>
<h1>Semantic terms</h1>
<ul>
    {% for term in bssFile.terms %}
        {% for key, value in term.items() %}
            <li>{{ key }}: {{ value }} </li>
        {% endfor %}
    {% endfor %}
</ul>
<h1>Datasets</h1>
        {% for dataset in bssFile.datasets %}
            <a href="{{ dataset.path }}">{{ dataset.path }}</a>
            {{ dataset.title }}
            <hr/>
        {% endfor %}
       </div> 
</main>

<script src="https://code.jquery.com/jquery-3.4.1.slim.min.js" integrity="sha384-J6qa4849blE2+poT4WnyKhv5vZF5SrPo0iEjwBvKU7imGFAV0wwj1yYfoRSJoZ+n" crossorigin="anonymous"></script>
<script src="https://cdn.jsdelivr.net/npm/popper.js@1.16.0/dist/umd/popper.min.js" integrity="sha384-Q6E9RHvbIyZFJoft+2mJbHaEWldlvI9IOYy5n3zV9zzTtmI3UksdQRVvoxMfooAo" crossorigin="anonymous"></script>
<script src="https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/js/bootstrap.min.js" integrity="sha384-wfSDF2E50Y2D1uUdj0O3uMBJnjuUD4Ih7YwaYd1iqfktj0Uod8GCExl3Og8ifwB6" crossorigin="anonymous"></script>
</body>
""")
        self.bssFile.add_entry("index.html", template.render(bssFile=self.bssFile))
        return ProcessorReturn.SUCCESSFUL
